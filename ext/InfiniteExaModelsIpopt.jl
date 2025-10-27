module InfiniteExaModelsIpopt

import InfiniteExaModels, NLPModelsIpopt, SolverCore, NLPModels

const _DefaultWallTime = 1.0e20
const _DefaultPrintLevel = 5
const _SilentPrintLevel = 0

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
    if backend.silent && get(prev_options, :print_level, _DefaultPrintLevel) != _SilentPrintLevel
        new_options[:print_level] = _SilentPrintLevel
    elseif !backend.silent && 
        get(prev_options, :print_level, _DefaultPrintLevel) == _SilentPrintLevel && 
        !haskey(options, :print_level)
        # If previously silent & not otherwise specified, restore to default
        new_options[:print_level] = _DefaultPrintLevel
    end
    # Process time limit setting
    if !isnan(backend.time_limit) && get(prev_options, :max_wall_time, NaN) != backend.time_limit
        new_options[:max_wall_time] = backend.time_limit
    elseif !haskey(options, :max_wall_time) && 
        isnan(backend.time_limit) && 
        get(prev_options, :max_wall_time, _DefaultWallTime) != _DefaultWallTime
        # If previously set & not otherwise specified, restore to default time limit
        new_options[:max_wall_time] = _DefaultWallTime
    end
    # Save updated options for more potential resolves
    merge!(prev_options, new_options)
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