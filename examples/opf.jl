using InfiniteExaModels
using InfiniteOpt, Distributions, NLPModelsIpopt, Ipopt, Random
using ExaModelsExamples

function main(filename = "pglib_opf_case14_ieee.m"; seed = 0, num_supports = 10)
    
    Random.seed!(seed)

    # Set the parameters
    n_θ = 3
    θ_nom = [0.; 60.; 10.]
    covar = [80. 0 0; 0 80. 0; 0 0 120.]
    ref = ExaModelsExamples.get_power_data_ref(filename)

    branch = [
        begin
            g,b = ExaModelsExamples.PowerModels.calc_branch_y(branch)
            tr,ti = ExaModelsExamples.PowerModels.calc_branch_t(branch)
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
    im = InfiniteModel(Ipopt.Optimizer)
    # set_silent(im)
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
        expect(
            sum(gen["cost"][1] * pg[i]^2 + gen["cost"][2] * pg[i] + gen["cost"][3] for (i, gen) in ref[:gen]),
            θ
        )
    )
    
    
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

        JuMP.@constraint(
            im,
            sum(p[a] for a in ref[:bus_arcs][i]) ==
            sum(pg[g] for g in ref[:bus_gens][i]) - sum(load["pd"] for load in bus_loads) -
            sum(shunt["gs"] for shunt in bus_shunts) * vm[i]^2
        )

        JuMP.@constraint(
            im,
            sum(q[a] for a in ref[:bus_arcs][i]) ==
            sum(qg[g] for g in ref[:bus_gens][i]) - sum(load["qd"] for load in bus_loads) +
            sum(shunt["bs"] for shunt in bus_shunts) * vm[i]^2
        )
    end

    # Create the ExaModel and solve both models to compare
    optimize!(im)
    
    @time em, mappings = exa_model(im)
    # result = ipopt(em)

    # # Get the answers
    # ed = [result.solution[mappings.finvar_mappings[v].i] for v in d]
    # id = value.(d)

    # # Print a report
    # println("\n--------------------------------------------")
    # println("               SUMMARY")
    # println("--------------------------------------------\n")
    # println("ExaModel Objective:      ", -result.objective) # change sign for maximization
    # println("InfiniteModel Objective: ", objective_value(im))
    # println("\nExaModel d:      ", ed)
    # println("InfiniteModel d: ", id)
end

main()
