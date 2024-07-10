using InfiniteExaModels
using InfiniteOpt, Distributions, Random
using PowerModels
using Downloads
using LinearAlgebra

function get_power_case(filename)
    if !isfile(filename)
        if !isdir("$(@__DIR__)/temp")
            mkdir("$(@__DIR__)/temp")
        end
        ff = joinpath(@__DIR__, "temp", filename)
        if !isfile(ff)
            @info "Downloading $filename"
            Downloads.download(
                "https://raw.githubusercontent.com/power-grid-lib/pglib-opf/dc6be4b2f85ca0e776952ec22cbd4c22396ea5a3/$filename",
                joinpath(@__DIR__, "temp", filename),
            )
            return joinpath(@__DIR__, "temp", filename)
        else
            return ff
        end
    else
        return filename
    end
end

function get_power_data_ref(filename)
    case = get_power_case(filename)
    data = PowerModels.parse_file(case)
    PowerModels.standardize_cost_terms!(data, order = 2)
    PowerModels.calc_thermal_limits!(data)
    return PowerModels.build_ref(data)[:it][:pm][:nw][0]
end

function opf(filename = "pglib_opf_case3_lmbd.m"; seed = 0, num_supports = 100, backend = nothing)
    
    ref = get_power_data_ref(filename)
    
    
    Random.seed!(seed)
    
    # Set the parameters
    nbus = length(ref[:bus])
    n_θ = nbus * 2
    θ_nom = zeros(n_θ)
    
    Pdcovar = [sum(load["pd"] for load in [ref[:load][l] for l in ref[:bus_loads][i]]) for (i, bus) in ref[:bus]]   
    Qdcovar = [sum(load["qd"] for load in [ref[:load][l] for l in ref[:bus_loads][i]]) for (i, bus) in ref[:bus]]   
    covar = (0.1 * [Pdcovar; Qdcovar]).^2
    
    branch = [
    begin
        g,b = PowerModels.calc_branch_y(branch)
        tr,ti = PowerModels.calc_branch_t(branch)
        f_bus = branch["f_bus"]
        t_bus = branch["t_bus"]
        (
        f_bus = f_bus,
        t_bus = t_bus,
        f_idx = (i, f_bus, t_bus),
        t_idx = (i, t_bus, f_bus),
        g = g,
        b = b,
        tr = tr,
        ti = ti,
        ttm = tr^2 + ti^2,
        g_fr = branch["g_fr"],
        b_fr = branch["b_fr"],
        g_to = branch["g_to"],
        b_to = branch["b_to"],
        angmin = branch["angmin"],
        angmax = branch["angmax"],
        rate_a = branch["rate_a"]
        )
    end
    for (i, branch) in ref[:branch]]
        
        # Define the infinite model
        im = InfiniteModel(backend)
        
        # first stage variables
        @variable(im, va0[i in keys(ref[:bus])])
        @variable(
            im,
            ref[:bus][i]["vmin"] <= vm0[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"],
            start = 1.0
        )
        @variable(
            im,
            ref[:gen][i]["pmin"] <= pg0[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"],
        )
        @variable(
            im,
            ref[:gen][i]["qmin"] <= qg0[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"],
        )
        
        @variable(
            im,
            -ref[:branch][l]["rate_a"] <=
            p0[(l, i, j) in ref[:arcs]] <=
            ref[:branch][l]["rate_a"]
        )
        @variable(
            im,
            -ref[:branch][l]["rate_a"] <=
            q0[(l, i, j) in ref[:arcs]] <=
            ref[:branch][l]["rate_a"]
        )
        
        # second stage variables
        @infinite_parameter(im, θ[i = 1:n_θ] ~ MvNormal(θ_nom, covar), num_supports = num_supports)
        
        @variable(im, va[i in keys(ref[:bus])], Infinite(θ))
        @variable(
            im,
            ref[:bus][i]["vmin"] <= vm[i in keys(ref[:bus])] <= ref[:bus][i]["vmax"], Infinite(θ),
            start = 1.0
        )
        @variable(
            im,
            ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"], Infinite(θ),
        )
        @variable(
            im,
            ref[:gen][i]["qmin"] <= qg[i in keys(ref[:gen])] <= ref[:gen][i]["qmax"], Infinite(θ),
        )
        
        @variable(
            im,
            -ref[:branch][l]["rate_a"] <=
            p[(l, i, j) in ref[:arcs]] <=
            ref[:branch][l]["rate_a"], Infinite(θ)
        )
        @variable(
            im,
            -ref[:branch][l]["rate_a"] <=
            q[(l, i, j) in ref[:arcs]] <=
            ref[:branch][l]["rate_a"], Infinite(θ)
        )
        
        @objective(
            im,
            Min,
            sum(gen["cost"][1] * pg0[i]^2 + gen["cost"][2] * pg0[i] + gen["cost"][3] for (i, gen) in ref[:gen])
        )
        
        # first stage constraints
        @constraint(im, [(i, bus) in ref[:ref_buses]], va0[i] == 0)
        @constraint(
            im,
            [br in branch],
            p0[br.f_idx] ==
            (br.g + br.g_fr) / br.ttm * vm0[br.f_bus]^2 +
            (-br.g * br.tr + br.b * br.ti) / br.ttm * (vm0[br.f_bus] * vm0[br.t_bus] * cos(va0[br.f_bus] - va0[br.t_bus])) +
            (-br.b * br.tr - br.g * br.ti) / br.ttm * (vm0[br.f_bus] * vm0[br.t_bus] * sin(va0[br.f_bus] - va0[br.t_bus]))
        )
        @constraint(
            im,
            [br in branch],
            q0[br.f_idx] ==
            -(br.b + br.b_fr) / br.ttm * vm0[br.f_bus]^2 -
            (-br.b * br.tr - br.g * br.ti) / br.ttm * (vm0[br.f_bus] * vm0[br.t_bus] * cos(va0[br.f_bus] - va0[br.t_bus])) +
            (-br.g * br.tr + br.b * br.ti) / br.ttm * (vm0[br.f_bus] * vm0[br.t_bus] * sin(va0[br.f_bus] - va0[br.t_bus]))
        )
        
        # To side of the branch flow
        @constraint(
            im,
            [br in branch],
            p0[br.t_idx] ==
            (br.g + br.g_to) * vm0[br.t_bus]^2 +
            (-br.g * br.tr - br.b * br.ti) / br.ttm * (vm0[br.t_bus] * vm0[br.f_bus] * cos(va0[br.t_bus] - va0[br.f_bus])) +
            (-br.b * br.tr + br.g * br.ti) / br.ttm * (vm0[br.t_bus] * vm0[br.f_bus] * sin(va0[br.t_bus] - va0[br.f_bus]))
        )
        @constraint(
            im,
            [br in branch],
            q0[br.t_idx] ==
            -(br.b + br.b_to) * vm0[br.t_bus]^2 -
            (-br.b * br.tr + br.g * br.ti) / br.ttm * (vm0[br.t_bus] * vm0[br.f_bus] * cos(va0[br.t_bus] - va0[br.f_bus])) +
            (-br.g * br.tr - br.b * br.ti) / br.ttm * (vm0[br.t_bus] * vm0[br.f_bus] * sin(va0[br.t_bus] - va0[br.f_bus]))
        )
        
        # Voltage angle difference limit
        @constraint(im, [br in branch], br.angmin <= va0[br.f_bus] - va0[br.t_bus] <= br.angmax)
        
        # Apparent power limit, from side and to side
        @constraint(im, [br in branch], p0[br.f_idx]^2 + q0[br.f_idx]^2 <= br.rate_a)
        @constraint(im, [br in branch], p0[br.t_idx]^2 + q0[br.t_idx]^2 <= br.rate_a)
        
        # Bottle-neck
        for (i, bus) in ref[:bus]
            bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
            bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
            
            @constraint(
                im,
                sum(p0[a] for a in ref[:bus_arcs][i]) ==
                sum(pg0[g] for g in ref[:bus_gens][i]) - sum(load["pd"] for load in bus_loads) -
                sum(shunt["gs"] for shunt in bus_shunts) * vm0[i]^2
            )
            
            @constraint(
                im,
                sum(q0[a] for a in ref[:bus_arcs][i]) ==
                sum(qg0[g] for g in ref[:bus_gens][i]) - sum(load["qd"] for load in bus_loads) +
                sum(shunt["bs"] for shunt in bus_shunts) * vm0[i]^2
            )
        end
        
        # second stage constraints
        @constraint(im, [(i, bus) in ref[:ref_buses]], va[i] == 0)
        @constraint(
            im,
            [br in branch],
            p[br.f_idx] ==
            (br.g + br.g_fr) / br.ttm * vm[br.f_bus]^2 +
            (-br.g * br.tr + br.b * br.ti) / br.ttm * (vm[br.f_bus] * vm[br.t_bus] * cos(va[br.f_bus] - va[br.t_bus])) +
            (-br.b * br.tr - br.g * br.ti) / br.ttm * (vm[br.f_bus] * vm[br.t_bus] * sin(va[br.f_bus] - va[br.t_bus]))
        )
        @constraint(
            im,
            [br in branch],
            q[br.f_idx] ==
            -(br.b + br.b_fr) / br.ttm * vm[br.f_bus]^2 -
            (-br.b * br.tr - br.g * br.ti) / br.ttm * (vm[br.f_bus] * vm[br.t_bus] * cos(va[br.f_bus] - va[br.t_bus])) +
            (-br.g * br.tr + br.b * br.ti) / br.ttm * (vm[br.f_bus] * vm[br.t_bus] * sin(va[br.f_bus] - va[br.t_bus]))
        )
        
        # To side of the branch flow
        @constraint(
            im,
            [br in branch],
            p[br.t_idx] ==
            (br.g + br.g_to) * vm[br.t_bus]^2 +
            (-br.g * br.tr - br.b * br.ti) / br.ttm * (vm[br.t_bus] * vm[br.f_bus] * cos(va[br.t_bus] - va[br.f_bus])) +
            (-br.b * br.tr + br.g * br.ti) / br.ttm * (vm[br.t_bus] * vm[br.f_bus] * sin(va[br.t_bus] - va[br.f_bus]))
        )
        @constraint(
            im,
            [br in branch],
            q[br.t_idx] ==
            -(br.b + br.b_to) * vm[br.t_bus]^2 -
            (-br.b * br.tr + br.g * br.ti) / br.ttm * (vm[br.t_bus] * vm[br.f_bus] * cos(va[br.t_bus] - va[br.f_bus])) +
            (-br.g * br.tr - br.b * br.ti) / br.ttm * (vm[br.t_bus] * vm[br.f_bus] * sin(va[br.t_bus] - va[br.f_bus]))
        )
        
        # Voltage angle difference limit
        @constraint(im, [br in branch], br.angmin <= va[br.f_bus] - va[br.t_bus] <= br.angmax)
        
        # Apparent power limit, from side and to side
        @constraint(im, [br in branch], p[br.f_idx]^2 + q[br.f_idx]^2 <= br.rate_a)
        @constraint(im, [br in branch], p[br.t_idx]^2 + q[br.t_idx]^2 <= br.rate_a)
        
        # Bottle-neck
        for (i, bus) in ref[:bus]
            bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
            bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]
            
            @constraint(
                im,
                sum(p[a] for a in ref[:bus_arcs][i]) ==
                θ[i] +
                sum(pg[g] for g in ref[:bus_gens][i]) -
                sum(load["pd"] for load in bus_loads) -
                sum(shunt["gs"] for shunt in bus_shunts) * vm[i]^2
            )
            
            @constraint(
                im,
                sum(q[a] for a in ref[:bus_arcs][i]) ==
                θ[nbus + i] +
                sum(qg[g] for g in ref[:bus_gens][i]) -
                sum(load["qd"] for load in bus_loads) +
                sum(shunt["bs"] for shunt in bus_shunts) * vm[i]^2
            )
        end
        
        # ramping constraints
        # TODO: these bounds are arbitrary. Need to check what's the practical value.
        @constraint(im, [(i, gen) in ref[:gen]],  - .1 * (gen["pmax"]-gen["pmin"]) <= pg0[i] - pg[i] <= .1 * (gen["pmax"]-gen["pmin"]))
        @constraint(im, [(i, gen) in ref[:gen]],  - .1 * (gen["qmax"]-gen["qmin"]) <= qg0[i] - qg[i] <= .1 * (gen["qmax"]-gen["qmin"]))
        
        return im
    end
    
    