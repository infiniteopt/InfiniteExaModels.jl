module InfiniteExaModelsIpopt

import InfiniteExaModels, NLPModelsIpopt, SolverCore

# Account for the silent and time limit settings
function _process_options(options, backend)
    prev_options = backend.prev_options
    # Get new or updated options to pass to solver
    new_options = Dict{Symbol, Any}(
        k => v 
        for (k, v) in options 
        if !haskey(prev_options, k) || prev_options[k] != v
    )
    # Process silent setting
    if backend.silent
        new_options[:print_level] = 0
    elseif iszero(get(prev_options, :print_level, 1)) && !haskey(new_options, :print_level)
        # If previously silent & not otherwise specified, restore to default
        new_options[:print_level] = 5
    end
    # Process time limit setting
    if !isnan(backend.time_limit) && get(prev_options, :max_wall_time, NaN) != backend.time_limit
        new_options[:max_wall_time] = backend.time_limit
    elseif !haskey(new_options, :max_wall_time) && isnan(backend.time_limit) && !isnan(get(prev_options, :max_wall_time, NaN))
        # If previously silent & not otherwise specified, restore to default time limit
        new_options[:max_wall_time] = 1.0e20
    end
    # Turn off warmstart
    if get(new_options, :warm_start_init_point, false) == "no"
        # Ensure no warmstarts are passed to solver
        delete!(new_options, :x0)
        delete!(new_options, :y0)
        delete!(new_options, :zL0)
        delete!(new_options, :zU0)
    end
    # Save updated options for more potential resolves
    backend.prev_options = new_options
    return new_options
end

# Setup the solver, solve it, and return the results
function InfiniteExaModels.initial_solve(
    type::Type{NLPModelsIpopt.IpoptSolver},
    backend,
    options
    )
    sol_options = _process_options(options, backend)
    backend.solver = type(backend.model)
    return SolverCore.solve!(backend.solver, backend.model; sol_options...)
end

# Prepare solver for resolve, solve it, and return the results
function InfiniteExaModels.resolve(
    solver::NLPModelsIpopt.IpoptSolver,
    backend,
    options
    )
    sol_options = _process_options(options, backend)
    SolverCore.reset!(solver, backend.model)
    if isnothing(backend.results)
        backend.results = SolverCore.GenericExecutionStats(backend.model)
    end
    return SolverCore.solve!(solver, backend.model, backend.results; sol_options...)
end

function InfiniteExaModels.warmstart_backend(
    backend::InfiniteExaModels.ExaTranscriptionBackend,
    solver::NLPModelsIpopt.IpoptSolver
    )
    results = backend.results
    options = backend.options
    # Add to options so that these can be passed onto the solver later on
    options[:x0] = results.solution
    options[:y0] = results.multipliers
    options[:zL0] = results.multipliers_L
    options[:zU0] = results.multipliers_U
    return
end
end