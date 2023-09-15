import InfiniteOpt, JuMP, ExaModels
import MathOptInterface as _MOI
import InfiniteOpt.TranscriptionOpt as _TO
using Base.Meta

struct MappingData
    # Mappings
    infvar_mappings::Dict{InfiniteOpt.GeneralVariableRef, ExaModels.Variable}
    finvar_mappings::Dict{InfiniteOpt.GeneralVariableRef, ExaModels.Var}
    constraint_mappings::Dict{InfiniteOpt.InfOptConstraintRef, ExaModels.Constraint}

    # Helpful metadata
    param_alias::Dict{InfiniteOpt.GeneralVariableRef, Symbol}
    group_alias::Vector{Symbol}
    base_itrs::Vector{Any}
    support_to_index::Dict{Tuple{Int, Union{Float64, Vector{Float64}}}, Int}
    
    # Default constructor
    function MappingData()
        return new(
            Dict{InfiniteOpt.GeneralVariableRef, ExaModels.Variable}(),
            Dict{InfiniteOpt.GeneralVariableRef, ExaModels.Var}(),
            Dict{InfiniteOpt.InfOptConstraintRef, ExaModels.Constraint}(),
            Dict{InfiniteOpt.GeneralVariableRef, Symbol}(),
            Symbol[],
            [],
            Dict{Tuple{Int, Union{Float64, Vector{Float64}}}, Int}(),
        )
    end
end

# Generate infinite parameter references based on basic indices
function _parameter_ref(
    model::InfiniteOpt.InfiniteModel,
    idx::InfiniteOpt.IndependentParameterIndex
    )
    return InfiniteOpt.GeneralVariableRef(model, idx.value, InfiniteOpt.IndependentParameterIndex)
end
function _parameter_ref(
    model::InfiniteOpt.InfiniteModel,
    idx::InfiniteOpt.DependentParametersIndex
    )
    num_params = length(model.dependent_params[idx].parameter_nums)
    return [InfiniteOpt.GeneralVariableRef(model, idx.value, InfiniteOpt.DependentParameterIndex, i) for i in 1:num_params]
end

# Create the support iterators for each infinite parameter group and add to the mapping data
function _build_base_iterators(
    data::MappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    # gather the individual infinite parameter groups
    param_indices = InfiniteOpt._param_object_indices(inf_model)
    param_groups = map(idx -> _parameter_ref(inf_model, idx), param_indices)
    # build the iterator for each parameter group
    for group in param_groups # group will either be singular ref or a vector of refs
        # generate all the symbols for created named tuples
        for pref in group
            data.param_alias[pref] = gensym()
        end
        aliases = map(pref -> data.param_alias[pref], group)
        itr_sym = gensym()
        push!(data.group_alias, itr_sym)
        # setup the supports (discretization points)
        InfiniteOpt.add_generative_supports(first(group))
        supps = keys(InfiniteOpt._parameter_supports(InfiniteOpt.dispatch_variable_ref(first(group))))
        group_idx = length(data.group_alias)
        for (i, s) in enumerate(supps)
            data.support_to_index[group_idx, s] = i
        end
        # create the iterator which is a tuple of named tuples
        if aliases isa Vector
            itr = Tuple((; itr_sym => i, zip(aliases, s)...) for (i, s) in enumerate(supps))
        else
            itr = Tuple((; itr_sym => i, aliases => s) for (i, s) in enumerate(supps))
        end
        # add the iterator to `data`
        push!(data.base_itrs, itr)
    end
    return
end

# Determine the bounds of an InfiniteOpt variable
function _get_variable_bounds(vref)
    lb = -Inf
    ub = Inf
    if JuMP.is_fixed(vref)
        lb = JuMP.fix_value(vref)
        ub = lb
    end
    if JuMP.has_lower_bound(vref)
        lb = JuMP.lower_bound(vref)
    end
    if JuMP.has_upper_bound(vref)
        ub = JuMP.upper_bound(vref)
    end
    return lb, ub
end

# All the finite variables from an InfiniteModel to a ExaCore
function _add_finite_variables(
    core::ExaModels.ExaCore, 
    data::MappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for vref in JuMP.all_variables(inf_model, InfiniteOpt.FiniteVariable)
        # check if is continuous
        if JuMP.is_binary(vref) || JuMP.is_integer(vref)
            error("Integer variables not supported by ExaModels.")
        end
        # get the start values and bounds
        start = JuMP.start_value(vref)
        start = isnothing(start) ? 0 : start
        lb, ub = _get_variable_bounds(vref)
        # make the ExaModel variable
        new_var = ExaModels.variable(core, 1, start = start, lvar = lb, uvar = ub)
        data.finvar_mappings[vref] = new_var[1]
    end
    return
end

# All the infinite variables from an InfiniteModel to a ExaCore
function _add_infinite_variables(
    core::ExaModels.ExaCore, 
    data::MappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for vref in JuMP.all_variables(inf_model, InfiniteOpt.InfiniteVariable) # TODO include derivative variables
        # check if is continuous
        if JuMP.is_binary(vref) || JuMP.is_integer(vref)
            error("Integer variables not supported by ExaModels.")
        end
        # retrieve the basic info
        var = InfiniteOpt._core_variable_object(vref)
        group_idxs = var.object_nums
        # iterate through the supports to get the start values 
        itrs = map(i -> data.base_itrs[i], group_idxs)
        dims = length.(itrs)
        starts = zeros(dims...)
        if !isnothing(InfiniteOpt.start_value_function(vref))
            for i in Iterators.product(itrs...)
                supp = [s for nt in i for s in Iterators.drop(values(nt), 2)] # TODO maybe store this for later 
                if var.is_vector_start
                    start = var.info.start(supp)
                else
                    start = var.info.start(Tuple(supp, var.parameter_refs)...)
                end
                starts[first.(i)...] = start
            end
        end
        # prepare the bounds
        lb, ub = _get_variable_bounds(vref)
        # create the ExaModels variable
        new_var = ExaModels.variable(core, dims...; start = starts, lvar = lb, uvar = ub)
        data.infvar_mappings[vref] = new_var
    end
    return 
end

# All the point variables from an InfiniteModel to a ExaCore
function _add_point_variables(
    core::ExaModels.ExaCore, 
    data::MappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for vref in JuMP.all_variables(inf_model, InfiniteOpt.PointVariable)
        # check if is continuous
        if JuMP.is_binary(vref) || JuMP.is_integer(vref)
            error("Integer variables not supported by ExaModels.")
        end
        # TODO account for changes in bounds and start
        # store the index mapping for the point variable
        ivref = InfiniteOpt.infinite_variable_ref(vref)
        raw_supp = InfiniteOpt.raw_parameter_values(vref)
        supp = Tuple(raw_supp, InfiniteOpt.raw_parameter_refs(ivref), use_indices = false)
        group_idxs = InfiniteOpt._object_numbers(ivref)
        idxs = Tuple(data.support_to_index[i, s] for (i, s) in zip(group_idxs, supp))
        data.finvar_mappings[vref] = data.infvar_mappings[ivref][idxs...]
    end
    return 
end

# Map variable references based on their underlying type
# Note: these all will have access input `p` which is named tuple from the iterator
function _map_variable(vref, data)
    return _map_variable(vref, vref.index_type, data)
end
function _map_variable(
    vref, 
    ::Type{V}, 
    data
    ) where V <: Union{InfiniteOpt.FiniteVariableIndex, InfiniteOpt.PointVariableIndex}
    return data.finvar_mappings[vref]
end
function _map_variable(
    vref, 
    ::Type{V}, 
    data
    ) where V <: InfiniteOpt.InfiniteVariableIndex
    group_idxs = InfiniteOpt._object_numbers(vref)
    ex = Expr(:ref, :($(data.infvar_mappings[vref])))
    append!(ex.args, :(p.$(data.group_alias[i])) for i in group_idxs)
    return ex
end
function _map_variable(
    vref, 
    ::Type{V}, 
    data
    ) where V <: InfiniteOpt.InfiniteParameterIndex
    return :(p.$(data.param_alias[vref]))
end
function _map_variable(
    vref, 
    ::Type{V}, 
    data
    ) where V <: InfiniteOpt.FiniteParameterIndex
    return InfiniteOpt.parameter_value(vref)
end
function _map_variable(
    vref, 
    IdxType, 
    data
    )
    error("Unable to add `$vref` to an ExaModel.")
end

# Construct a Julia function that is compatible with a generator for building ExaModels
function _make_expr_gen_func(expr, data)
    expr_code = InfiniteOpt.map_expression_to_ast(v -> _map_variable(v, data), expr) # TODO add operator mapper argument
    func_code = :(p -> $expr_code) # p is a namedtuple with all the parameters
    f = eval(func_code)
    return x -> Base.invokelatest(f, x) # needed to avoid world age problems
end

# Extract the constraint bounds from an MOI set
function _get_constr_bounds(set::_MOI.LessThan)
    return -Inf, set.upper
end
function _get_constr_bounds(set::_MOI.GreaterThan)
    return set.lower, Inf
end
function _get_constr_bounds(::_MOI.EqualTo)
    return 0, 0
end
function _get_constr_bounds(set::_MOI.Interval)
    return set.lower, set.upper
end
function _get_constr_bounds(set)
    error("Constraint set `$set` is not compatible with ExaModels.")
end

# Add all the constraints from an InfiniteModel to an ExaCore
function _add_constraints(
    core::ExaModels.ExaCore, 
    data::MappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for cref in JuMP.all_constraints(inf_model)
        # skip if the constraint is a variable bound or type
        InfiniteOpt._is_info_constraint(cref) && continue
        # parse the basic information
        constr = JuMP.constraint_object(cref)
        expr = JuMP.jump_function(constr)
        set = JuMP.moi_set(constr)
        group_idxs = InfiniteOpt._object_numbers(cref)
        # make the expression generator
        func = _make_expr_gen_func(expr, data)
        if isempty(group_idxs) # we have a finite constraint
            itr = ((;),)
        elseif length(group_idxs) == 1 # we only depend on one independent infinite parameter
            itr = data.base_itrs[first(group_idxs)]
        else # we depend on multiple independent infinite parameters
            itrs = map(i -> data.base_itrs[i], group_idxs)
            itr = Tuple(merge(i...) for i in Iterators.product(itrs...))
        end
        # TODO account for DomainRestrictions
        gen = Base.Generator(func, itr)
        # get the constraint bounds
        lb, ub = _get_constr_bounds(set)
        # create the ExaModels constraint
        con = ExaModels.constraint(core, gen, lcon = lb, ucon = ub)
        data.constraint_mappings[cref] = con
    end
    return
end

# Fill an ExaCore and MappingData with an InfiniteModel
function build_exa_core!(
    core::ExaModels.ExaCore, 
    data::MappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    _build_base_iterators(data, inf_model)
    InfiniteOpt.expand_all_measures!(inf_model)
    _add_finite_variables(core, data, inf_model)
    _add_infinite_variables(core, data, inf_model)
    _add_point_variables(core, data, inf_model)
    _add_constraints(core, data, inf_model)
    return
end

# Create an ExaModel from an InfiniteModel
function exa_model(inf_model::InfiniteOpt.InfiniteModel)
    core = ExaModels.ExaCore()
    data = MappingData()
    build_exa_core!(core, data, inf_model)
    return ExaModels.ExaModel(core), data 
end
