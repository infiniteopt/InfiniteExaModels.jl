using ExaModels, DelimitedFiles, MadNLPGPU, CUDA
include("quadrotor.jl")
include("opf.jl")
include("pandemic.jl")
include("utils.jl")

function transopt_madnlp_test(im, filename)
    set_optimizer_attribute(im, "output_file", joinpath(@__DIR__, "logs","$(filename)_madnlp.log"))
    total_time = @elapsed optimize!(im)
    obj = objective_value(im)
    status = primal_status(im)
    nvar, ncon, sol_time, ad_time = madnlp_stats(joinpath(@__DIR__, "logs","$(filename)_madnlp.log"))
    return nvar, ncon, obj, status, total_time, sol_time, ad_time
end

function infexa_madnlp_test(im, filename; backend = nothing)
    total_time = @elapsed begin 
        em, _ = exa_model(im; backend = backend)
        result = madnlp(em; output_file = joinpath(@__DIR__, "logs","$(filename)_madnlp.log"))
    end
    obj = result.objective
    status = result.status
    nvar, ncon, sol_time, ad_time = madnlp_stats(joinpath(@__DIR__, "logs","$(filename)_madnlp.log"))
    return nvar, ncon, obj, status, total_time, sol_time, ad_time
end

function run_cases(im_func, kwarg_dict_list; prerun = true, backend = nothing)
    # add directories if needed
    if !isdir("$(@__DIR__)/results")
        mkdir("$(@__DIR__)/results")
    end
    if !isdir("$(@__DIR__)/logs")
        mkdir("$(@__DIR__)/logs")
    end
    # run everything once for the JIT compiler
    opt_dict = (
        "ExaModelsMOI" => () -> ExaModels.MadNLPOptimizer(backend),
        "InfiniteExaModels" => nothing
        )
    first_kwargs = first(kwarg_dict_list)
    if prerun
        for (opt_name, opt) in opt_dict
            im = im_func(; opt = opt, first_kwargs...)
            if opt_name == "InfiniteExaModels"
                infexa_madnlp_test(im, "compile_run", backend = backend)
            else
                transopt_madnlp_test(im, "compile_run")
            end
        end
    end

    # initialize a table to save the results to
    kwarg_strs = sort!(string.(keys(first_kwargs)))
    opts = []
    nvars = []
    ncons = []
    objs = []
    statuses = []
    total_times = []
    sol_times = []
    ad_times = []
    for (opt_name, opt) in opt_dict
        for d in kwarg_dict_list
            im = im_func(; opt = opt, d...)
            filename = string(opt_name, "_", join([d[Symbol(k)] for k in kwarg_strs], "_"))
            if opt_name == "InfiniteExaModels"
                output = infexa_madnlp_test(im, filename, backend = backend)
            else
                output = transopt_madnlp_test(im, filename)
            end
            nvar, ncon, obj, status, total_time, sol_time, ad_time = output
            push!(opts, opt_name)
            push!(nvars, nvar)
            push!(ncons, ncon)
            push!(objs, obj)
            push!(statuses, status)
            push!(total_times, total_time)
            push!(sol_times, sol_time)
            push!(ad_times, ad_time)
        end
    end

    # make a table to save the results to
    nrows = length(kwarg_dict_list) * length(opt_dict) + 1
    ncols = length(first_kwargs) + 8
    table = Matrix{String}(undef, nrows, ncols)
    table[1, :] = append!(copy(kwarg_strs), ["framework", "nvar", "ncon", "objective", "status", "total_time", "solve_time", "ad_time"])
    kwarg_matrix = hcat([repeat([d[Symbol(k)] for d in kwarg_dict_list], length(opt_dict)) for k in kwarg_strs]...)
    table[2:end, :] = string.(hcat(kwarg_matrix, opts, nvars, ncons, objs, statuses, total_times, sol_times, ad_times))

    # save the matrix as a CSV
    open("$(@__DIR__)/results/$(im_func)_madnlp_results.csv", "w") do io
        writedlm(io, table, ",")
    end
end

# Run the opf study
num_supports_list = [1000, 2000, 4000, 8000, 16000]
settings = [Dict(:num_supports => n) for n in num_supports_list]
run_cases(opf, settings, backend = CUDABackend())

# Run the quadcopter study
num_supports_list = [1000, 2000, 4000, 8000, 16000]
settings = [Dict(:num_supports => n) for n in num_supports_list]
run_cases(quad, settings, backend = CUDABackend())

# Run the pandemic study
num_supports_list = [(25, 4), (50, 4), (100, 4), (100, 8), (100, 128)]
settings = [Dict(:num_supports => nt, :num_scenarios => nxi) for (nt, nxi) in num_supports_list]
run_cases(pandemic, settings, backend = CUDABackend())