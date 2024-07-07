"""

"""
struct ExaMappingData
    # Mappings
    infvar_mappings::Dict{InfiniteOpt.GeneralVariableRef, ExaModels.Variable}
    finvar_mappings::Dict{InfiniteOpt.GeneralVariableRef, ExaModels.Var}
    constraint_mappings::Dict{InfiniteOpt.InfOptConstraintRef, ExaModels.Constraint}

    # Helpful metadata
    param_alias::Dict{InfiniteOpt.GeneralVariableRef, Symbol}
    group_alias::Vector{Symbol}
    base_itrs::Vector{Any}
    support_to_index::Dict{Tuple{Int, Union{Float64, Vector{Float64}}}, Int}
    semivar_info::Dict{InfiniteOpt.GeneralVariableRef, Tuple{ExaModels.Variable, Vector{Any}}}
    pfunc_info::Dict{InfiniteOpt.GeneralVariableRef, Tuple{Symbol, <:Array}}
    
    # Default constructor
    function ExaMappingData()
        return new(
            Dict{InfiniteOpt.GeneralVariableRef, ExaModels.Variable}(),
            Dict{InfiniteOpt.GeneralVariableRef, ExaModels.Var}(),
            Dict{InfiniteOpt.InfOptConstraintRef, ExaModels.Constraint}(),
            Dict{InfiniteOpt.GeneralVariableRef, Symbol}(),
            Symbol[],
            [],
            Dict{Tuple{Int, Union{Float64, Vector{Float64}}}, Int}(),
            Dict{InfiniteOpt.GeneralVariableRef, Tuple{ExaModels.Variable, Vector{Any}}}(),
            Dict{InfiniteOpt.GeneralVariableRef, Tuple{Symbol, <:Array}}(),
        )
    end
end

# Store default options that differ from standard JSO defaults
const _DefaultOptions = Dict{Symbol, Any}(:verbose => 1, :max_time => Inf)

"""

"""
mutable struct ExaTranscriptionBackend{B} <: InfiniteOpt.AbstractTransformationBackend
    model::Union{Nothing, ExaModels.ExaModel}
    backend::B
    solver::Union{Nothing, SolverCore.AbstractOptimizationSolver}
    options::Dict{Symbol, Any}
    results::Union{Nothing, SolverCore.AbstractExecutionStats}
    data::ExaMappingData
end

"""

"""
function ExaTranscriptionBackend(; backend = nothing)
    return ExaTranscriptionBackend(
        nothing,
        backend,
        nothing,
        copy(_DefaultOptions),
        nothing,
        ExaMappingData()
        )
end
function ExaTranscriptionBackend(
    solver_type::Union{Type{<:SolverCore.AbstractOptimizationSolver}, MOI.OptimizerWithAttributes};
    kwargs...
    )
    backend = ExaTranscriptionBackend(; kwargs...)
    JuMP.set_optimizer(backend, solver_type)
    return backend
end
function ExaTranscriptionBackend(solver_type; kwargs...)
    error("`ExaTranscriptionBackend`s only support solver constructors of type " *
          " `SolverCore.AbstractOptimizationSolver` (e.g., `NLPModelsIpopt.IpoptSolver`).")
end

# Extend Base.empty!
function Base.empty!(backend::ExaTranscriptionBackend)
    backend.model = nothing
    backend.solver = nothing
    backend.results = nothing
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
    ) # TODO maybe add kwargs
    empty!(backend)
    backend.model = ExaModels.ExaModel(model, backend.data; backend = backend.backend)
    return
end

## Solver settings
# Symbol
function JuMP.get_attribute(backend::ExaTranscriptionBackend, attr::Symbol)
    haskey(backend.options, attr) || error("Attribute `$attr` not found.")
    return backend.options[attr]
end
function JuMP.set_attribute(backend::ExaTranscriptionBackend, attr::Symbol, value)
    return backend.options[attr] = value
end

# String
function JuMP.get_attribute(backend::ExaTranscriptionBackend, attr::AbstractString)
    return JuMP.get_attribute(backend, Symbol(attr))
end
function JuMP.set_attribute(backend::ExaTranscriptionBackend, attr::AbstractString, value)
    return JuMP.set_attribute(backend, Symbol(attr), value)
end

# MOI.RawOptimizerAttribute
function JuMP.get_attribute(backend::ExaTranscriptionBackend, attr::MOI.RawOptimizerAttribute)
    return JuMP.get_attribute(backend, attr.name)
end
function JuMP.set_attribute(backend::ExaTranscriptionBackend, attr::MOI.RawOptimizerAttribute, value)
    return JuMP.set_attribute(backend, attr.name, value)
end

# MOI.Silent
function JuMP.get_attribute(backend::ExaTranscriptionBackend, ::MOI.Silent)
    return JuMP.get_attribute(backend, :verbose) == 0
end
function JuMP.set_attribute(backend::ExaTranscriptionBackend, ::MOI.Silent, value::Bool)
    return JuMP.set_attribute(backend, :verbose, Int(!value))
end

# MOI.TimeLimitSec
function JuMP.get_attribute(backend::ExaTranscriptionBackend, ::MOI.TimeLimitSec)
    return JuMP.get_attribute(backend, :max_time)
end
function JuMP.set_attribute(backend::ExaTranscriptionBackend, ::MOI.TimeLimitSec, value::Real)
    return JuMP.set_attribute(backend, :max_time, Float64(value))
end
function JuMP.set_attribute(backend::ExaTranscriptionBackend, ::MOI.TimeLimitSec, ::Nothing)
    return JuMP.set_attribute(backend, :max_time, _DefaultOptions[:max_time])
end

# MOI.SolverName
function JuMP.get_attribute(backend::ExaTranscriptionBackend, ::MOI.SolverName)
    return string(get(backend.options, :solver, "No solver attached"))
end

# Optimizer specification
function JuMP.set_optimizer(
    backend::ExaTranscriptionBackend,
    solver_type::Type{<:SolverCore.AbstractOptimizationSolver}
    )
    # clear out all previous solver specific settings and set the solver
    merge!(empty!(backend.options), _DefaultOptions)
    backend.options[:solver] = solver_type
    # clear out the old solver and results
    backend.solver = nothing
    backend.results = nothing
    return
end
function JuMP.set_optimizer(
    backend::ExaTranscriptionBackend,
    solver_type::MOI.OptimizerWithAttributes
    )
    JuMP.set_optimizer(backend, solver_type.optimizer_constructor)
    for (attr, value) in solver_type.params
        JuMP.set_attribute(backend, attr, value)
    end
    return
end

# Fallback for translating the default option nomenclature to the solver
function translate_option(p::Pair, solver_type)
    return p
end

# Extend Optimize!
function JuMP.optimize!(backend::ExaTranscriptionBackend)
    if isnothing(backend.solver)
        haskey(backend.options, :solver) || throw(JuMP.NoOptimizer())
        solver_type = backend.options[:solver]
        options = (translate_option(p, solver_type) for p in backend.options if p[1] != :solver)
        backend.solver = solver_type(backend.model; options...)
    else
        SolverCore.reset!(backend.solver, backend.model)
    end
    backend.results = SolverCore.solve!(backend.solver, backend.model)
    return
end

# InfiniteOpt mapping methods
function InfiniteOpt.transformation_variable(
    vref::InfiniteOpt.GeneralVariableRef,
    backend::ExaTranscriptionBackend
    ) # TODO add kwargs
    data = backend.data
    haskey(data.infvar_mappings, vref) && return data.infvar_mappings[vref]
    haskey(data.finvar_mappings, vref) && return data.finvar_mappings[vref]
    error("A mapping for `$vref` in the transformation backend not found.")
end

## Result queries
# MOI.ResultCount
function JuMP.get_attribute(backend::ExaTranscriptionBackend, ::MOI.ResultCount)
    return isnothing(backend.results) ? 0 : 1
end
