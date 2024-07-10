module InfiniteExaModelsMadNLP

import InfiniteExaModels, MadNLP
import MathOptInterface as _MOI

function InfiniteExaModels.translate_option(p::Pair, ::Type{MadNLP.MadNLPSolver})
    if p[1] == :verbose
        return :print_level => isone(p[2]) ? 5 : 0
    else
        return p
    end
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