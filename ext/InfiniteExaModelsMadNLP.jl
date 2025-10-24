module InfiniteExaModelsMadNLP

import InfiniteExaModels, MadNLP, SolverCore, NLPModels
import MathOptInterface as _MOI

const _DefaultWallTime = 1.0E6
const _DefaultPrintLevel = MadNLP.INFO
const _SilentPrintLevel = MadNLP.ERROR

# Set new solver options + update existing ones
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
    type::Type{MadNLP.MadNLPSolver},
    backend,
    options
    )
    sol_options = _process_options(options, backend)
    backend.solver = type(backend.model; sol_options...)
    return MadNLP.solve!(backend.solver)
end

# Prepare solver for resolve, solve it, and return the results
function InfiniteExaModels.resolve(
    solver::MadNLP.MadNLPSolver,
    backend,
    options
    )
    sol_options = _process_options(options, backend)
    # Update solver logger print level (for resolves)
    if haskey(sol_options, :print_level)
        backend.solver.logger.print_level = sol_options[:print_level]
    end
    return MadNLP.solve!(backend.model, solver; sol_options...)
end

# Standard JSO statuses to MOI.TerminationStatusCode
const _TerminationMappings = Dict{MadNLP.Status, _MOI.TerminationStatusCode}(
    MadNLP.SOLVE_SUCCEEDED => _MOI.LOCALLY_SOLVED,
    MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL => _MOI.ALMOST_LOCALLY_SOLVED,
    MadNLP.SEARCH_DIRECTION_BECOMES_TOO_SMALL => _MOI.SLOW_PROGRESS,
    MadNLP.DIVERGING_ITERATES => _MOI.INFEASIBLE_OR_UNBOUNDED,
    MadNLP.INFEASIBLE_PROBLEM_DETECTED => _MOI.LOCALLY_INFEASIBLE,
    MadNLP.MAXIMUM_ITERATIONS_EXCEEDED => _MOI.ITERATION_LIMIT,
    MadNLP.MAXIMUM_WALLTIME_EXCEEDED => _MOI.TIME_LIMIT,
    MadNLP.INITIAL => _MOI.OPTIMIZE_NOT_CALLED,
    MadNLP.RESTORATION_FAILED => _MOI.NUMERICAL_ERROR,
    MadNLP.INVALID_NUMBER_DETECTED => _MOI.INVALID_MODEL,
    MadNLP.ERROR_IN_STEP_COMPUTATION => _MOI.NUMERICAL_ERROR,
    MadNLP.NOT_ENOUGH_DEGREES_OF_FREEDOM => _MOI.INVALID_MODEL,
    MadNLP.USER_REQUESTED_STOP => _MOI.INTERRUPTED,
    MadNLP.INTERNAL_ERROR => _MOI.OTHER_ERROR,
    MadNLP.INVALID_NUMBER_OBJECTIVE => _MOI.INVALID_MODEL,
    MadNLP.INVALID_NUMBER_GRADIENT => _MOI.INVALID_MODEL,
    MadNLP.INVALID_NUMBER_CONSTRAINTS => _MOI.INVALID_MODEL,
    MadNLP.INVALID_NUMBER_JACOBIAN => _MOI.INVALID_MODEL,
    MadNLP.INVALID_NUMBER_HESSIAN_LAGRANGIAN => _MOI.INVALID_MODEL,
)

# Standard JSO statuses to MOI.ResultStatusCode
const _ResultMappings = Dict{MadNLP.Status, _MOI.ResultStatusCode}(
    MadNLP.SOLVE_SUCCEEDED => _MOI.FEASIBLE_POINT,
    MadNLP.SOLVED_TO_ACCEPTABLE_LEVEL => _MOI.NEARLY_FEASIBLE_POINT,
    MadNLP.INFEASIBLE_PROBLEM_DETECTED => _MOI.INFEASIBLE_POINT
)

# Default translation with standard JSO codes
function InfiniteExaModels.translate_termination_status(::MadNLP.MadNLPSolver, status)
    return get(_TerminationMappings, status, _MOI.OTHER_ERROR)
end

# Default translation with standard JSO codes
function InfiniteExaModels.translate_result_status(::MadNLP.MadNLPSolver, status)
    return get(_ResultMappings, status, _MOI.UNKNOWN_RESULT_STATUS)
end

end