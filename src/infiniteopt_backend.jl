"""

"""
struct ExaMappingData
    # Mappings
    infvar_mappings::Dict{InfiniteOpt.GeneralVariableRef, ExaModels.Variable}
    finvar_mappings::Dict{InfiniteOpt.GeneralVariableRef, ExaModels.Var}
    constraint_mappings::Dict{InfiniteOpt.InfOptConstraintRef, ExaModels.Constraint}
    param_mappings::Dict{InfiniteOpt.GeneralVariableRef, ExaModels.Parameter}

    # Helpful metadata
    param_alias::Dict{InfiniteOpt.GeneralVariableRef, Symbol}
    group_alias::Vector{Symbol}
    base_itrs::Vector{Any}
    support_labels::Vector{Vector{Set{DataType}}}
    has_internal_supps::Vector{Bool}
    support_to_index::Dict{Tuple{Int, Union{Float64, Vector{Float64}}}, Int}
    # Semi-infinite variable info
    semivar_info::Dict{
        InfiniteOpt.GeneralVariableRef,
        Tuple{
            Union{ExaModels.Variable, ExaModels.Parameter},
            Vector{Any}
        }
    }
    
    # Default constructor
    function ExaMappingData()
        return new(
            Dict{InfiniteOpt.GeneralVariableRef, ExaModels.Variable}(),
            Dict{InfiniteOpt.GeneralVariableRef, ExaModels.Var}(),
            Dict{InfiniteOpt.InfOptConstraintRef, ExaModels.Constraint}(),
            Dict{InfiniteOpt.GeneralVariableRef, ExaModels.Parameter}(),
            Dict{InfiniteOpt.GeneralVariableRef, Symbol}(),
            Symbol[],
            [],
            Vector{Set{DataType}}[],
            Bool[],
            Dict{Tuple{Int, Union{Float64, Vector{Float64}}}, Int}(),
            Dict{
                InfiniteOpt.GeneralVariableRef,
                Tuple{
                    Union{ExaModels.Variable, ExaModels.Parameter},
                    Vector{Any}
                }
            }(),
        )
    end
end

"""

"""
mutable struct ExaTranscriptionBackend{B} <: InfiniteOpt.AbstractTransformationBackend
    core::Union{Nothing, ExaModels.ExaCore}
    model::Union{Nothing, ExaModels.ExaModel}
    backend::B
    solver::Any
    prev_options::Dict{Symbol, Any}
    options::Dict{Symbol, Any}
    silent::Bool
    time_limit::Float64
    results::Union{Nothing, SolverCore.AbstractExecutionStats}
    solve_time::Float64
    data::ExaMappingData
end

"""

"""
function ExaTranscriptionBackend(; backend = nothing)
    return ExaTranscriptionBackend(
        nothing,
        nothing,
        backend,
        nothing,
        Dict{Symbol, Any}(),
        Dict{Symbol, Any}(),
        false,
        NaN,
        nothing,
        NaN,
        ExaMappingData()
        )
end
function ExaTranscriptionBackend(solver_type; kwargs...)
    backend = ExaTranscriptionBackend(; kwargs...)
    JuMP.set_optimizer(backend, solver_type)
    return backend
end

# Extend Base.empty! (leave the options alone)
function Base.empty!(backend::ExaTranscriptionBackend)
    backend.model = nothing
    backend.solver = nothing
    backend.results = nothing
    backend.solve_time = NaN
    backend.data = ExaMappingData()
    return backend
end

# Basic InfiniteOpt transformation backend extensions
InfiniteOpt.transformation_model(backend::ExaTranscriptionBackend) = backend.model
InfiniteOpt.transformation_data(backend::ExaTranscriptionBackend) = backend.data

# Build out the backend
function InfiniteOpt.build_transformation_backend!(
    model::InfiniteOpt.InfiniteModel,
    backend::ExaTranscriptionBackend
    )
    empty!(backend)
    backend.core = ExaModels.ExaCore(model, backend.data; backend = backend.backend)
    backend.model = ExaModels.ExaModel(backend.core)
end

## Solver settings
# Symbol
function JuMP.get_attribute(backend::ExaTranscriptionBackend, attr::Symbol)
    haskey(backend.options, attr) || error("Attribute `$attr` not found.")
    return backend.options[attr]
end
function JuMP.set_attribute(backend::ExaTranscriptionBackend, attr::Symbol, value)
    backend.solve_time = NaN
    backend.options[attr] = value
    return
end

# MOI.RawOptimizerAttribute
function JuMP.get_attribute(
    backend::ExaTranscriptionBackend,
    attr::_MOI.RawOptimizerAttribute
    )
    return JuMP.get_attribute(backend, Symbol(attr.name))
end
function JuMP.set_attribute(
    backend::ExaTranscriptionBackend,
    attr::_MOI.RawOptimizerAttribute,
    value
    )
    return JuMP.set_attribute(backend, Symbol(attr.name), value)
end

# MOI.Silent
function JuMP.get_attribute(
    backend::ExaTranscriptionBackend,
    ::_MOI.Silent
    )
    return backend.silent
end
function JuMP.set_attribute(
    backend::ExaTranscriptionBackend,
    ::_MOI.Silent,
    value::Bool
    )
    backend.silent = value
    return
end

# MOI.TimeLimitSec
function JuMP.get_attribute(
    backend::ExaTranscriptionBackend,
    ::_MOI.TimeLimitSec
    )
    limit = backend.time_limit
    return isnan(limit) ? nothing : limit
end
function JuMP.set_attribute(
    backend::ExaTranscriptionBackend,
    ::_MOI.TimeLimitSec,
    value::Real
    )
    backend.time_limit = Float64(value)
    return
end
function JuMP.set_attribute(
    backend::ExaTranscriptionBackend,
    ::_MOI.TimeLimitSec,
    ::Nothing
    )
    backend.time_limit = NaN
    return
end

# MOI.SolverName
function JuMP.get_attribute(backend::ExaTranscriptionBackend, ::_MOI.SolverName)
    return string(get(backend.options, :solver, "No solver attached"))
end

# Optimizer specification
function JuMP.set_optimizer(
    backend::ExaTranscriptionBackend,
    solver_type
    )
    # clear out all previous solver specific settings and set the solver
    empty!(backend.options)
    JuMP.set_attribute(backend, :solver, solver_type)
    backend.solver = nothing
    return
end
function JuMP.set_optimizer(
    backend::ExaTranscriptionBackend,
    solver_type::_MOI.OptimizerWithAttributes
    )
    JuMP.set_optimizer(backend, solver_type.optimizer_constructor)
    for (attr, value) in solver_type.params
        JuMP.set_attribute(backend, attr, value)
    end
    return
end

# TODO add JSO fallbacks
function initial_solve end
function resolve end

# Extend Optimize!
function JuMP.optimize!(backend::ExaTranscriptionBackend)
    haskey(backend.options, :solver) || throw(JuMP.NoOptimizer())
    solver_type = backend.options[:solver]
    options = filter(p -> p[1] != :solver, backend.options)
    backend.solve_time = @elapsed begin
        backend.results = if isnothing(backend.solver)
            initial_solve(solver_type, backend, options)
        else
            resolve(backend.solver, backend, options)
        end
    end
    return backend.results
end

# Check whether a mapping exists
function _check_mapping(vref::InfiniteOpt.GeneralVariableRef, backend)
    data = backend.data
    if !haskey(data.infvar_mappings, vref) && !haskey(data.finvar_mappings, vref) && !haskey(data.param_mappings, vref)
        error("A mapping for `$vref` in the transformation backend not found.")
    end
    return
end
function _check_mapping(cref::InfiniteOpt.InfOptConstraintRef, backend)
    if !haskey(backend.data.constraint_mappings, cref)
        error("A mapping for `$cref` in the transformation backend not found.")
    end
    return
end

# InfiniteOpt object mapping methods
function InfiniteOpt.transformation_variable(
    vref::InfiniteOpt.GeneralVariableRef,
    backend::ExaTranscriptionBackend
    ) # TODO add keywords (e.g., label)
    _check_mapping(vref, backend)
    data = backend.data
    haskey(data.infvar_mappings, vref) && return data.infvar_mappings[vref]
    haskey(data.finvar_mappings, vref) && return data.finvar_mappings[vref]
    return data.param_mappings[vref]
end
# TODO find a way to support expressions
function InfiniteOpt.transformation_constraint(
    cref::InfiniteOpt.InfOptConstraintRef,
    backend::ExaTranscriptionBackend
    ) # TODO add keywords (e.g., label)
    _check_mapping(cref, backend)
    return backend.data.constraint_mappings[cref]
end

# Extract the supports for a given variable/constraint reference
function _get_supports(ref, pref_vt, backend)
    group_idxs = InfiniteOpt.parameter_group_int_indices(ref)
    itrs = map(i -> backend.data.base_itrs[i], group_idxs)
    dims = length.(itrs)
    SuppType = typeof(Tuple(zeros(length(pref_vt)), pref_vt))
    supps = Array{SuppType, length(dims)}(undef, dims...)
    for i in Iterators.product(itrs...)
        supp = [s for nt in i for s in Iterators.drop(values(nt), 1)]
        supps[first.(i)...] = Tuple(supp, pref_vt)
    end
    return supps
end

# Helper function to filter output array based on pref label
# TODO this will need to be changed onced DomainRestrictions are added
function _label_filter(arr, ref, label, data)
    # check any filtering is needed
    label == InfiniteOpt.All && return arr
    group_idxs = InfiniteOpt.parameter_group_int_indices(ref)
    if label == InfiniteOpt.PublicLabel && 
        !any(data.has_internal_supps[group_idxs])
        return arr
    end
    # filter the array based on the desired label
    idxs = (map(s -> any(l -> l <: label, s), sets) for sets in data.support_labels)
    return arr[idxs...]
end

# InfiniteOpt support mapping methods
function InfiniteOpt.variable_supports(
    vref::Union{
        InfiniteOpt.InfiniteVariableRef,
        InfiniteOpt.SemiInfiniteVariableRef, 
        InfiniteOpt.DerivativeRef,
        InfiniteOpt.ParameterFunctionRef,
        InfiniteOpt.MeasureRef
        },
    backend::ExaTranscriptionBackend;
    label::DataType = InfiniteOpt.PublicLabel
    )
    pref_vt = InfiniteOpt.raw_parameter_refs(vref)
    supps = _get_supports(vref, pref_vt, backend)
    return _label_filter(supps, vref, label, backend.data)
end
# TODO find a way to support expressions
function InfiniteOpt.constraint_supports(
    cref::InfiniteOpt.InfOptConstraintRef,
    backend::ExaTranscriptionBackend;
    label::DataType = InfiniteOpt.PublicLabel
    )
    prefs = InfiniteOpt.parameter_refs(cref)
    pref_vt = InfiniteOpt.Collections.VectorTuple(prefs)
    supps = _get_supports(cref, pref_vt, backend)
    return _label_filter(supps, cref, label, backend.data)
end

# Check that a result is available
function _check_results_available(backend)
    if isnothing(backend.results)
        error("No solution available to query.")
    end
    return
end

# Standard JSO statuses to MOI.TerminationStatusCode
const _TerminationMappings = Dict{Symbol, _MOI.TerminationStatusCode}(
    :first_order => _MOI.LOCALLY_SOLVED,
    :acceptable => _MOI.ALMOST_LOCALLY_SOLVED,
    :small_step => _MOI.SLOW_PROGRESS,
    :infeasible => _MOI.INFEASIBLE_OR_UNBOUNDED,
    :unbounded => _MOI.INFEASIBLE_OR_UNBOUNDED,
    :max_iter => _MOI.ITERATION_LIMIT,
    :max_time => _MOI.TIME_LIMIT,
    :user => _MOI.INTERRUPTED,
    :exception => _MOI.OTHER_ERROR,
    :stalled => _MOI.OTHER_ERROR,
    :max_eval => _MOI.OTHER_LIMIT,
    :neg_pred => _MOI.OTHER_ERROR,
    :not_desc => _MOI.OTHER_ERROR,
)

# Standard JSO statuses to MOI.ResultStatusCode
const _ResultMappings = Dict{Symbol, _MOI.ResultStatusCode}(
    :first_order => _MOI.FEASIBLE_POINT,
    :acceptable => _MOI.NEARLY_FEASIBLE_POINT,
    :infeasible => _MOI.INFEASIBLE_POINT,
)

# Default translation with standard JSO codes
function translate_termination_status(solver, status)
    return get(_TerminationMappings, status, _MOI.OTHER_ERROR)
end

# Default translation with standard JSO codes
function translate_result_status(solver, status)
    return get(_ResultMappings, status, _MOI.UNKNOWN_RESULT_STATUS)
end

## Result queries
# MOI.ResultCount
function JuMP.get_attribute(
    backend::ExaTranscriptionBackend,
    ::_MOI.ResultCount
    )
    return isnothing(backend.results) ? 0 : 1
end

# MOI.RawStatusString
function JuMP.get_attribute(
    backend::ExaTranscriptionBackend,
    ::_MOI.RawStatusString
    )
    isnothing(backend.results) && return "optimize not called"
    return string(backend.results.status)
end

# MOI.TerminationStatus
function JuMP.get_attribute(
    backend::ExaTranscriptionBackend,
    ::_MOI.TerminationStatus
    )
    isnothing(backend.results) && return _MOI.OPTIMIZE_NOT_CALLED
    return translate_termination_status(backend.solver, backend.results.status)
end

# MOI.PrimalStatus/MOI.DualStatus
function JuMP.get_attribute(
    backend::ExaTranscriptionBackend,
    ::Union{_MOI.PrimalStatus, _MOI.DualStatus}
    )
    isnothing(backend.results) && return _MOI.NO_SOLUTION
    return translate_result_status(backend.solver, backend.results.status)
end

# MOI.SolveTimeSec
function JuMP.get_attribute(
    backend::ExaTranscriptionBackend,
    ::_MOI.SolveTimeSec
    )
    _check_results_available(backend)
    return backend.solve_time
end

# MOI.ObjectiveValue
function JuMP.get_attribute(
    backend::ExaTranscriptionBackend,
    ::_MOI.ObjectiveValue
    )
    _check_results_available(backend)
    return backend.results.objective
end

# map_value
function InfiniteOpt.map_value(
    vref::InfiniteOpt.GeneralVariableRef,
    backend::ExaTranscriptionBackend;
    label::DataType = InfiniteOpt.PublicLabel
    )
    _check_results_available(backend)
    var = InfiniteOpt.transformation_variable(vref, backend)
    vals = ExaModels.solution(backend.results, var)
    return _label_filter(vals, vref, label, backend.data)
end

# TODO find a way to support expressions/constraints

# Process variable bounds to give the right dual value
_get_domain_dual(mL, mU, ::_MOI.LessThan) = min.(mL .- mU, 0.0)
_get_domain_dual(mL, mU, ::_MOI.GreaterThan) = max.(mL .- mU, 0.0)
_get_domain_dual(mL, mU, ::_MOI.EqualTo) = mL .- mU

# map_dual
function InfiniteOpt.map_dual(
    cref::InfiniteOpt.InfOptConstraintRef,
    backend::ExaTranscriptionBackend;
    label::DataType = InfiniteOpt.PublicLabel
    )
    _check_results_available(backend)
    if InfiniteOpt.is_variable_domain_constraint(cref)
        con = JuMP.constraint_object(cref)
        var = InfiniteOpt.transformation_variable(JuMP.jump_function(con))
        set = JuMP.moi_set(con)
        mL = ExaModels.multipliers_L(backend.results, var) 
        mU = ExaModels.multipliers_U(backend.results, var)
        duals = _get_domain_dual(mL, mU, set)
    else
        con = InfiniteOpt.transformation_constraint(cref, backend)
        duals = -1.0 * ExaModels.multipliers(backend.results, con)
    end
    return _label_filter(duals, cref, label, backend.data)
end

# For updating parameters
function InfiniteOpt.update_parameter_value(
    backend::ExaTranscriptionBackend,
    pref::InfiniteOpt.FiniteParameterRef,
    value::Real
    )
    data = backend.data
    core = backend.core
    pref = InfiniteOpt.GeneralVariableRef(pref)
    # Check if the mapping exists
    haskey(data.param_mappings, pref) || return false
    # Update the value in the ExaCore, which updates the ExaModel too
    ExaModels.set_parameter!(
        core,
        data.param_mappings[pref],
        [value]
    )
    return true
end

# For updating parameter functions
function InfiniteOpt.update_parameter_value(
    backend::ExaTranscriptionBackend,
    pfref::InfiniteOpt.ParameterFunctionRef,
    value::Function
    )
    data = backend.data
    core = backend.core
    supps = InfiniteOpt.variable_supports(pfref, backend; label = InfiniteOpt.All)
    pfref = InfiniteOpt.GeneralVariableRef(pfref)
    # Check if the mapping exists
    haskey(data.param_mappings, pfref) || return false
    # Generate the new parameter function values
    vals = map(supp -> value(supp...), supps)
    # Update the values in the ExaCore, which updates the ExaModel too
    param = data.param_mappings[pfref]
    ExaModels.set_parameter!(core, param, vals)
    return true
end

# TODO add map_shadow_price, map_optimizer_index, map_reduced_cost
