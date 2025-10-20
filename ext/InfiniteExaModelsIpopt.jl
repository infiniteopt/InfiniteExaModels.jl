module InfiniteExaModelsIpopt

import InfiniteExaModels, NLPModelsIpopt, SolverCore

# Account for the silent and time limit settings
function _process_options(options, backend)
    prev_options = backend.prev_options

    # Process silent setting
    if backend.silent
        options[:print_level] = 0
    elseif !haskey(options, :print_level)
        if haskey(prev_options, :print_level) && iszero(prev_options[:print_level])
            # If previously silent, set to default
            options[:print_level] = 5
        end
    end

    # Process time limit setting
    if !isnan(backend.time_limit)
        options[:max_wall_time] = backend.time_limit
    else
        # Delete time limit if previously set
        delete!(options, :max_wall_time)
    end

    # Get new or updated options to pass to solver
    new_options = Dict(
        k => v 
        for (k, v) in options 
        if !haskey(prev_options, k) || prev_options[k] != v
    )
    # Save updated options for more potential resolves
    backend.prev_options = copy(new_options)
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
    return SolverCore.solve!(solver, backend.model, backend.results; sol_options...)
end

end