module InfiniteExaModelsIpopt

import InfiniteExaModels, NLPModelsIpopt, SolverCore

# Account for the silent and time limit settings
function _process_options(options, backend)
    if isnothing(backend.prev_options)
        # Process silent setting
        if backend.silent
            options[:print_level] = 0
        elseif !haskey(options, :print_level)
            # Default print level
            options[:print_level] = 5
        end
        # Process time limit setting
        if !isnan(backend.time_limit)
            options[:max_wall_time] = backend.time_limit
        end

        # Save options for potential resolve
        backend.prev_options = options
    else
        prev = backend.prev_options
        # Process silent setting in options
        if backend.silent
            # Updating to silent
            options[:print_level] = 0
        elseif !haskey(options, :print_level) && prev[:print_level] == 0
            # Update to default if print level not specified
            options[:print_level] = 5
        end
        # Process time limit setting in options
        if !isnan(backend.time_limit)
            options[:max_wall_time] = backend.time_limit
        else
            options[:max_wall_time] = NaN
        end

        # Keep only new or updated options to pass into solve!
        prev = backend.prev_options
        is_new = pair -> !haskey(prev, pair.first) || prev[pair.first] != pair.second
        options = filter(is_new, options)

        # Save updated options for more potential resolves
        merge!(backend.prev_options, options)
    end
    return
end

# Setup the solver, solve it, and return the results
function InfiniteExaModels.initial_solve(
    type::Type{NLPModelsIpopt.IpoptSolver},
    backend,
    options
    )
    _process_options(options, backend)
    backend.solver = type(backend.model)
    return SolverCore.solve!(backend.solver, backend.model; options...)
end

# Prepare solver for resolve, solve it, and return the results
function InfiniteExaModels.resolve(
    solver::NLPModelsIpopt.IpoptSolver,
    backend,
    options
    )
    _process_options(options, backend)
    # TODO handle removal of silence/time settings
    SolverCore.reset!(solver, backend.model)
    return SolverCore.solve!(solver, backend.model; options...)
end

end