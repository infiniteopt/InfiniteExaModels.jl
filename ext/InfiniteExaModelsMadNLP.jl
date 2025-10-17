module InfiniteExaModelsMadNLP

import InfiniteExaModels, MadNLP, SolverCore
import MathOptInterface as _MOI

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
        options[:print_level] = MadNLP.ERROR
    end
    if !isnan(backend.time_limit)
        options[:max_wall_time] = backend.time_limit
    end
    return
end

# Setup the solver, solve it, and return the results
function InfiniteExaModels.initial_solve(
    type::Type{MadNLP.MadNLPSolver},
    backend,
    options
    )
    _process_options(options, backend)
    backend.solver = type(backend.model; options...)
    return MadNLP.solve!(backend.solver)
end

# Prepare solver for resolve, solve it, and return the results
function InfiniteExaModels.resolve(
    solver::MadNLP.MadNLPSolver,
    backend,
    options
    )
    _process_options(options, backend)
    # TODO handle removal of silence/time settings
    # TODO fix this to handle options correctly (might require fix with MadNLP)
    return MadNLP.solve!(backend.model, solver; options...)
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