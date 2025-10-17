module InfiniteExaModelsIpopt

import InfiniteExaModels, NLPModelsIpopt, SolverCore

# Account for the silent and time limit settings
function _process_options(options, backend)
        if isnothing(backend.prev_options)
        # Save the options for potential resolve
        backend.prev_options = options
    else
        # Keep only new or updated options to pass into solve!
        prev = backend.prev_options
        is_new = pair -> !haskey(prev, pair.first) || prev[pair.first] != pair.second
        options = filter(is_new, options)

        # Save updated options for more potential resolves
        merge!(backend.prev_options, options)
    end
    if backend.silent
        options[:print_level] = 0
    end
    if !isnan(backend.time_limit)
        options[:max_wall_time] = backend.time_limit
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