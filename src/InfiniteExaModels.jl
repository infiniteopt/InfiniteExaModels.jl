module InfiniteExaModels

import InfiniteOpt, JuMP, ExaModels
import MathOptInterface as _MOI
import InfiniteOpt.TranscriptionOpt as _TO

export MappingData, build_exa_core!, exa_model

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
    semivar_info::Dict{InfiniteOpt.GeneralVariableRef, Tuple{ExaModels.Variable, Vector{Any}}}
    pfunc_info::Dict{InfiniteOpt.GeneralVariableRef, Tuple{Symbol, <:Array}}
    
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
            Dict{InfiniteOpt.GeneralVariableRef, Tuple{ExaModels.Variable, Vector{Any}}}(),
            Dict{InfiniteOpt.GeneralVariableRef, Tuple{Symbol, <:Array}}(),
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
            data.param_alias[pref] = if group isa Vector
                Symbol("dp$(pref.raw_index)$(pref.param_index)")
            else
                Symbol("ip$(pref.raw_index)")
            end
        end
        aliases = map(pref -> data.param_alias[pref], group)
        itr_sym = Symbol("group_idx$(length(data.group_alias)+1)")
        push!(data.group_alias, itr_sym)
        # setup the supports (discretization points)
        InfiniteOpt.add_generative_supports(first(group))
        supps = keys(InfiniteOpt._parameter_supports(InfiniteOpt.dispatch_variable_ref(first(group))))
        group_idx = length(data.group_alias)
        for (i, s) in enumerate(supps)
            data.support_to_index[group_idx, s] = i
        end
        # create the iterator which is a vector of named tuples
        itr = [(; itr_sym => i, zip(aliases, s)...) for (i, s) in enumerate(supps)]
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

# Add all the finite variables from an InfiniteModel to a ExaCore
function _add_finite_variables(
    core::ExaModels.ExaCore, 
    data::MappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for vref in JuMP.all_variables(inf_model, InfiniteOpt.FiniteVariable)
        # check if is continuous
        if JuMP.is_binary(vref) || JuMP.is_integer(vref)
            error("Integer variables are not supported by ExaModels.")
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

# Add all the infinite variables from an InfiniteModel to a ExaCore
function _add_infinite_variables(
    core::ExaModels.ExaCore, 
    data::MappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    ivrefs = JuMP.all_variables(inf_model, InfiniteOpt.InfiniteVariable)
    drefs = InfiniteOpt.all_derivatives(inf_model)
    for vref in append!(ivrefs, drefs)
        # check if is continuous
        if JuMP.is_binary(vref) || JuMP.is_integer(vref)
            error("Integer variables are not supported by ExaModels.")
        end
        # retrieve the basic info
        var = InfiniteOpt._core_variable_object(vref)
        group_idxs = InfiniteOpt._object_numbers(vref)
        # iterate through the supports to get the start values 
        itrs = map(i -> data.base_itrs[i], group_idxs)
        dims = length.(itrs)
        starts = zeros(dims...)
        if !isnothing(InfiniteOpt.start_value_function(vref))
            for i in Iterators.product(itrs...)
                supp = [s for nt in i for s in Iterators.drop(values(nt), 1)]
                if var.is_vector_start
                    start = var.info.start(supp)
                else
                    prefs = InfiniteOpt.raw_parameter_refs(vref)
                    start = var.info.start(Tuple(supp, prefs)...)
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

# Process all the parameter function from an InfiniteModel and prepare for use in a ExaCore
function _add_parameter_functions(
    data::MappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )  
    for (idx, obj) in inf_model.param_functions
        # gather the basic information
        pfref = InfiniteOpt._make_variable_ref(inf_model, idx)
        group_idxs = obj.func.object_nums
        prefs = obj.func.parameter_refs # a VectorTuple
        # compute the value for each support combination and store
        itrs = map(i -> data.base_itrs[i], group_idxs)
        dims = length.(itrs)
        vals = zeros(dims...)
        for i in Iterators.product(itrs...)
            supp = [s for nt in i for s in Iterators.drop(values(nt), 1)]
            vals[first.(i)...] = obj.func.func(Tuple(supp, prefs)...)
        end
        alias = Symbol("pf$(idx.value)")
        data.pfunc_info[pfref] = (alias, vals)
    end
    return
end

# Helper function for processing semi-infinite variables
function _process_semi_infinite_var(vref, data)
    # get basic info from InfiniteOpt
    ivref = InfiniteOpt.infinite_variable_ref(vref)
    orig_group_idxs = InfiniteOpt._object_numbers(ivref)
    raw_prefs = InfiniteOpt.raw_parameter_refs(ivref) # type `InfiniteOpt.VectorTuple`
    group_idxs = InfiniteOpt._object_numbers(vref)
    eval_supps = InfiniteOpt.eval_supports(vref)
    # create metadata vector `indexing` 
    indexing = Vector{Any}(undef, length(orig_group_idxs))
    for (i, g) in enumerate(orig_group_idxs)
        if g in group_idxs
            indexing[i] = data.group_alias[g] # get the group alias for indexing
        elseif isnothing(raw_prefs.indices[i])
            pref_ind = first(raw_prefs.ranges[i])
            supp = eval_supps[pref_ind]
            indexing[i] = data.support_to_index[g, supp] # store the support index
        else
            pref_inds = raw_prefs.ranges[i]
            supp = [eval_supps[i] for i in pref_inds]
            indexing[i] = data.support_to_index[g, supp] # store the support index
        end
    end
    # store the desired information
    if ivref.index_type == InfiniteOpt.ParameterFunctionIndex
        indexing[isa.(indexing, Symbol)] .= Colon()
        alias = Symbol("rpf$(vref.raw_index)")
        data.pfunc_info[vref] = (alias, data.pfunc_info[ivref][2][indexing...])
    else
        data.semivar_info[vref] = (data.infvar_mappings[ivref], indexing)
    end
    return
end

# Add all the semi-infinite variables from an InfiniteModel to a ExaCore
# In other words, create helpful metadata to be used by `_map_variable`
function _add_semi_infinite_variables(
    data::MappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for vref in JuMP.all_variables(inf_model, InfiniteOpt.SemiInfiniteVariable)
        # check if is continuous
        # if JuMP.is_binary(vref) || JuMP.is_integer(vref)
        #     error("Integer variables are not supported by ExaModels.")
        # end
        # TODO account for changes in bounds and start
        # collect the basic information fields from InfiniteOpt and save
        _process_semi_infinite_var(vref, data)
    end
    return 
end

# Helper function for processing point variables
function _process_point_var(vref, data)
    ivref = InfiniteOpt.infinite_variable_ref(vref)
    raw_supp = InfiniteOpt.raw_parameter_values(vref)
    supp = Tuple(raw_supp, InfiniteOpt.raw_parameter_refs(ivref), use_indices = false)
    group_idxs = InfiniteOpt._object_numbers(ivref)
    idxs = Tuple(data.support_to_index[i, s] for (i, s) in zip(group_idxs, supp))
    return data.infvar_mappings[ivref][idxs...]
end

# Add all the point variables from an InfiniteModel to a ExaCore
function _add_point_variables(
    data::MappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for vref in JuMP.all_variables(inf_model, InfiniteOpt.PointVariable)
        # check if is continuous
        if JuMP.is_binary(vref) || JuMP.is_integer(vref)
            error("Integer variables are not supported by ExaModels.")
        end
        # TODO account for changes in bounds and start
        # store the index mapping for the point variable
        data.finvar_mappings[vref] = _process_point_var(vref, data)
    end
    return 
end

# Map variable references based on their underlying type
function _map_variable(
    vref, 
    ::Type{InfiniteOpt.FiniteVariableIndex}, 
    itr_par,
    data
    )
    return data.finvar_mappings[vref]
end
function _map_variable(
    vref, 
    ::Type{InfiniteOpt.PointVariableIndex}, 
    itr_par,
    data
    )
    if haskey(data.finvar_mappings, vref) 
        return data.finvar_mappings[vref]
    else
        var = _process_point_var(vref, data)
        data.finvar_mappings[vref] = var
        return var
    end
end
function _map_variable(
    vref, 
    ::Type{V}, 
    itr_par,
    data
    ) where V <: Union{InfiniteOpt.InfiniteVariableIndex, InfiniteOpt.DerivativeIndex}
    group_idxs = InfiniteOpt._object_numbers(vref)
    idx_pars = (itr_par[data.group_alias[i]] for i in group_idxs)
    return data.infvar_mappings[vref][idx_pars...]
end
function _map_variable(
    vref, 
    ::Type{InfiniteOpt.SemiInfiniteVariableIndex}, 
    itr_par,
    data
    )
    if !haskey(data.semivar_info, vref) && !haskey(data.pfunc_info, vref)
        _process_semi_infinite_var(vref, data)
    end
    if haskey(data.semivar_info, vref)
        ivar, inds = data.semivar_info[vref]
        idx_pars = (i isa Int ? i : itr_par[i] for i in inds)
        return ivar[idx_pars...]
    else # we have a reduced parameter function
        return itr_par[data.pfunc_info[vref][1]]
    end
end
function _map_variable(
    vref, 
    ::Type{<:InfiniteOpt.InfiniteParameterIndex}, 
    itr_par,
    data
    )
    return itr_par[data.param_alias[vref]]
end
function _map_variable(
    vref, 
    ::Type{InfiniteOpt.FiniteParameterIndex}, 
    itr_par,
    data
    )
    return InfiniteOpt.parameter_value(vref)
end
function _map_variable(
    vref, 
    ::Type{InfiniteOpt.ParameterFunctionIndex}, 
    itr_par,
    data
    )
    return itr_par[data.pfunc_info[vref][1]]
end
function _map_variable(
    vref, 
    IdxType, 
    itr_par,
    data
    )
    error("Unable to add `$vref` to an ExaModel, it's index type `$IdxType`" *
          " is not yet supported by InfiniteExaModels.")
end

# Store the mappings between JuMP's operator names and the functions supported by ExaModels
const _op_mappings = Dict{Symbol, Function}(
    :+ => +,
    :- => -,
    :* => *,
    :/ => /,
    :^ => ^,
    :inv => inv,
    :sqrt => sqrt,
    :cbrt => cbrt,
    :abs => abs,
    :abs2 => abs2,
    :exp => exp,
    :exp2 => exp2,
    :log => log,
    :log2 => log2,
    :log10 => log10,
    :log1p => log1p,
    :sin => sin,
    :cos => cos,
    :tan => tan,
    :asin => asin,
    :acos => acos,
    :csc => csc,
    :sec => sec,
    :cot => cot,
    :atan => atan,
    :acot => acot,
    :sind => sind,
    :cosd => cosd,
    :tand => tand,
    :cscd => cscd,
    :secd => secd,
    :cotd => cotd,
    :atand => atand,
    :acotd => acotd,
    :sinh => sinh,
    :cosh => cosh,
    :tanh => tanh,
    :csch => csc,
    :sech => sech,
    :coth => coth,
    :atanh => atanh,
    :acoth => acoth,
    # TODO add the remaining JuMP operators
)

# Map a nonlinear function symbol to the underlying function
function _nl_op(s::Symbol)
    if !haskey(_op_mappings, s)
        error("`InfiniteExaModel`s does not support the nonlinear operator `$s`.")
    end
    return _op_mappings[s]
end

# Convert as InfiniteOpt expression into a ExaModel expression using the ParIndexed `itr_par`
# TODO ensure that this never returns a raw constant at the very end, maybe use ExaModels.Null...
function _exafy(vref::InfiniteOpt.GeneralVariableRef, itr_par, data)
    return _map_variable(vref, vref.index_type, itr_par, data)
end
function _exafy(c::Number, itr_par, data)
    return c
end
function _exafy(
    aff::JuMP.GenericAffExpr{C, InfiniteOpt.GeneralVariableRef}, 
    itr_par, 
    data
    ) where {C}
    c = JuMP.constant(aff)
    if !isempty(aff.terms)
        ex = sum(begin
            v_ex = _exafy(v, itr_par, data) 
            isone(c) ? v_ex : c * v_ex
            end for (c, v) in JuMP.linear_terms(aff)
            )
        return iszero(c) ? ex : ex + c
    else
        return c
    end
end
function _exafy(
    quad::JuMP.GenericQuadExpr{C, InfiniteOpt.GeneralVariableRef}, 
    itr_par, 
    data
    ) where {C}
    aff = _exafy(quad.aff, itr_par, data)
    if !isempty(quad.terms)
        ex = sum(begin 
            if v1 == v2
                v_ex = _exafy(v1, itr_par, data) 
                isone(c) ? abs2(v_ex) : c * abs2(v_ex)
            else
                v1_ex = _exafy(v1, itr_par, data) 
                v2_ex = _exafy(v2, itr_par, data) 
                isone(c) ? v1_ex * v2_ex : c * v1_ex * v2_ex
            end
            end for (c, v1, v2) in JuMP.quad_terms(quad)
            )
        return iszero(quad.aff) ? ex : ex + aff
    else
        return aff
    end
end
function _exafy(
    nl::JuMP.GenericNonlinearExpr{InfiniteOpt.GeneralVariableRef}, 
    itr_par, 
    data
    )
    return _nl_op(nl.head)((_exafy(a, itr_par, data) for a in nl.args)...)
end

# Extract the constraint bounds from an MOI set
function _get_constr_bounds(set::_MOI.LessThan)
    return -Inf, set.upper
end
function _get_constr_bounds(set::_MOI.GreaterThan)
    return set.lower, Inf
end
function _get_constr_bounds(set::_MOI.EqualTo)
    return set.value, set.value
end
function _get_constr_bounds(set::_MOI.Interval)
    return set.lower, set.upper
end
function _get_constr_bounds(set)
    error("Constraint set `$set` is not compatible with InfiniteExaModels.")
end

# Helper function to augment an expression function and iterator to incorporate parameter functions
# NOTE: _make_expr_ast must be called on `expr` before this is called
function _add_parameter_functions_to_itr(expr, itr, data)
    vrefs = InfiniteOpt._all_function_variables(expr)
    filter!(vrefs) do v
        if v.index_type == InfiniteOpt.ParameterFunctionIndex
            return true
        elseif v.index_type == InfiniteOpt.SemiInfiniteVariableIndex
            return InfiniteOpt.infinite_variable_ref(v).index_type == InfiniteOpt.ParameterFunctionIndex
        else
            return false
        end
    end
    # TODO figure out how to only update the iterator once
    for pfref in vrefs
        group_idxs = InfiniteOpt._object_numbers(pfref)
        aliases = [data.group_alias[g] for g in group_idxs]
        alias, vals = data.pfunc_info[pfref]
        itr = [(; p..., alias => vals[values(p[aliases])...]) for p in itr]
    end
    return itr
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
        if isempty(inf_model.constraints[JuMP.index(cref)].measure_indices)
            expr = JuMP.jump_function(constr)
        else
            @warn "Constrained measures can lead to poor performance with ExaModels."
            expr = InfiniteOpt.expand_measures(JuMP.jump_function(constr), inf_model)
        end
        set = JuMP.moi_set(constr)
        group_idxs = InfiniteOpt._object_numbers(cref)
        # prepare the iterator of NamedTuples (contains support values, iterator values, and constants from parameter functions)
        if isempty(group_idxs) # we have a finite constraint
            itr = [(;)]
        elseif length(group_idxs) == 1 # we only depend on one independent infinite parameter
            itr = data.base_itrs[first(group_idxs)]
        else # we depend on multiple independent infinite parameters
            itrs = map(i -> data.base_itrs[i], group_idxs)
            itr = [merge(i...) for i in Iterators.product(itrs...)]
        end
        # TODO account for DomainRestrictions
        if InfiniteOpt.has_domain_restrictions(cref)
            error("`DomainRestrictions` are not currently supported by InfiniteExaModels.")
        end
        # update the iterator to include parameter function values
        if !isempty(data.pfunc_info)
            itr = _add_parameter_functions_to_itr(expr, itr, data)
        end
        # create the ExaModels expression tree based on expr
        itr_par = ExaModels.Par(typeof(first(itr)))
        em_expr = _exafy(expr, itr_par, data)
        # get the constraint bounds
        lb, ub = _get_constr_bounds(set)
        # create the ExaModels constraint
        con = ExaModels.constraint(core, em_expr, itr, lcon = lb, ucon = ub)
        data.constraint_mappings[cref] = con
    end
    return
end

# Create the parameter iterator for finite difference methods
function _finite_diff_itr(
    method::InfiniteOpt.FiniteDifference{InfiniteOpt.Backward}, 
    srt_itr,
    p_alias
    )
    offset = zip(Iterators.drop(srt_itr, 1), srt_itr) # current, previous
    itr = [merge(p1, (; dt = p1[p_alias] - p2[p_alias], i1 = p1[1], i2 = p2[1])) for (p1, p2) in offset]
    return method.add_boundary_constraint ? itr : itr[1:end-1]
end
function _finite_diff_itr(
    method::InfiniteOpt.FiniteDifference{InfiniteOpt.Forward}, 
    srt_itr,
    p_alias
    )
    offset = zip(Iterators.drop(srt_itr, 1), srt_itr) # next, current
    itr = [merge(p2, (; dt = p1[p_alias] - p2[p_alias], i1 = p1[1], i2 = p2[1])) for (p1, p2) in offset]
    return method.add_boundary_constraint ? itr : itr[2:end]
end
function _finite_diff_itr(
    ::InfiniteOpt.FiniteDifference{InfiniteOpt.Central}, 
    srt_itr,
    p_alias
    )
    offset = zip(Iterators.drop(srt_itr, 2), Iterators.drop(srt_itr, 1), srt_itr) # next, current, prev
    return [merge(pc, (; dt = p1[p_alias] - p2[p_alias], i1 = p1[1], i2 = p2[1])) for (p1, pc, p2) in offset]
end
function _finite_diff_itr(method, srt_itr, p_alias)
    error("Unsupported finite difference approximation technique `$(method.technique)`.")
end

# Add auxiliary equations to make deriviative variables well-defined
function _add_deriv_equations(core, dref, data, method::InfiniteOpt.FiniteDifference)
    # gather the variables
    vref = InfiniteOpt.derivative_argument(dref) 
    pref = InfiniteOpt.operator_parameter(dref)
    # gather the needed infinite parameter data
    group_idxs = InfiniteOpt._object_numbers(vref)
    pref_group = InfiniteOpt._object_number(pref)
    # prepare the iterator for the equations
    p_alias = data.param_alias[pref]
    base_itr = data.base_itrs[pref_group]
    if length(first(base_itr)) == 2 # we have a independent parameter (the supports are sorted already)
        srt_itr = base_itr
    else # we have a dependent parameter (supports are not sorted)
        srt_itr = sort(base_itr, by = p -> p[p_alias])
    end
    pref_itr = _finite_diff_itr(method, srt_itr, p_alias)
    itrs = [g == pref_group ? pref_itr : data.base_itrs[g] for g in group_idxs]
    itr = [merge(i...) for i in Iterators.product(itrs...)]
    # make the expression function for the generator
    dvar = data.infvar_mappings[dref]
    d_inds = map(g -> data.group_alias[g], group_idxs)
    if vref.index_type == InfiniteOpt.SemiInfiniteVariableIndex
        ivar, inds = data.semivar_info[vref]
        g_alias = data.group_alias[pref_group]
        inds1 = Tuple(i == g_alias ? :i1 : i for i in inds)
        inds2 = Tuple(i == g_alias ? :i2 : i for i in inds)
        f = p -> dvar[values(p[d_inds])...] * p.dt - ivar[(i isa Int ? i : p[i] for i in inds1)...] + ivar[(i isa Int ? i : p[i] for i in inds2)...]
    elseif vref.index_type == InfiniteOpt.InfiniteVariableIndex
        ivar = data.infvar_mappings[vref]
        inds1 = Tuple(g == pref_group ? :i1 : data.group_alias[g] for g in group_idxs)
        inds2 = Tuple(g == pref_group ? :i2 : data.group_alias[g] for g in group_idxs)
        f = p -> dvar[values(p[d_inds])...] * p.dt - ivar[values(p[inds1])...] + ivar[values(p[inds2])...]
    else # TODO try to make measures work...
        error("Derivatives that act on references with index type `$(vref.index_type)` are not supported by InfiniteExaModels.")
    end
    # add the equations to core
    gen = Base.Generator(f, itr)
    ExaModels.constraint(core, gen) # TODO add mapping?
    return 
end
function _add_deriv_equations(
    core, 
    dref, 
    data, 
    method::InfiniteOpt.OrthogonalCollocation{InfiniteOpt.MeasureToolbox.GaussLobatto}
    )
    # gather the variables
    vref = InfiniteOpt.derivative_argument(dref) 
    pref = InfiniteOpt.operator_parameter(dref)
    # gather the needed infinite parameter data
    group_idxs = InfiniteOpt._object_numbers(vref)
    pref_group = InfiniteOpt._object_number(pref)
    pref_alias = data.group_alias[pref_group]
    aliases = Tuple(data.group_alias[g] for g in group_idxs)
    other_idxs = filter(i -> i != pref_group, group_idxs)
    # prepare the iterator for base equations
    n_nodes = method.num_nodes - 1 # we want the # of nodes in each element (excluding the lb)
    n_supps = length(data.base_itrs[pref_group])
    pref_itr = [(i1 = n, i2 = lb) for (n, lb) in zip(2:n_supps, repeat(1:n_nodes:n_supps-1, inner = n_nodes))]
    other_itrs = Tuple(data.base_itrs[g] for g in other_idxs)
    base_itr = vec([merge(i...) for i in Iterators.product(other_itrs..., pref_itr)])
    # create the base constraints (the summation portion excluded)
    if vref.index_type == InfiniteOpt.SemiInfiniteVariableIndex
        ivar, inds = data.semivar_info[vref]
        inds1 = Tuple(i == pref_alias ? :i1 : i for i in inds)
        inds2 = Tuple(i == pref_alias ? :i2 : i for i in inds)
        f = p -> -ivar[(i isa Int ? i : p[i] for i in inds1)...] + ivar[(i isa Int ? i : p[i] for i in inds2)...]
    elseif vref.index_type == InfiniteOpt.InfiniteVariableIndex
        ivar = data.infvar_mappings[vref]
        inds1 = Tuple(a == pref_alias ? :i1 : a for a in aliases)
        inds2 = Tuple(a == pref_alias ? :i2 : a for a in aliases)
        f = p -> -ivar[values(p[inds1])...] + ivar[values(p[inds2])...]
    else # TODO try to make measures work...
        error("Derivatives that act on references with index type `$(vref.index_type)` are not supported by InfiniteExaModels.")
    end
    base_con = ExaModels.constraint(core, Base.Generator(f, base_itr))
    # augment the base constraints with the summation terms
    all_supps = InfiniteOpt.supports(pref, label = InfiniteOpt.All)
    num_repeats = reduce(*, length.(other_itrs), init = 1)
    linear_idxs = LinearIndices(Tuple(1:length(data.base_itrs[g]) for g in group_idxs))
    itr_len = (n_supps-1) * n_nodes * num_repeats
    aug_itr = Vector{NamedTuple{(:i, :c, :idx), Tuple{Int, Float64, Int}}}(undef, itr_len)
    counter = 1
    for i in 1:n_nodes:n_supps-1
        # collect the supports
        lb = all_supps[i]
        i_nodes = all_supps[i+1:i+n_nodes] .- lb
        # build the matrices in memory order (column-wise using the transpose form)
        M1t = Matrix{Float64}(undef, n_nodes, n_nodes)
        M2t = Matrix{Float64}(undef, n_nodes, n_nodes)
        for j in eachindex(i_nodes) # support index
            for k in eachindex(i_nodes) # polynomial index
                M1t[k, j] = k * i_nodes[j]^(k-1)
                M2t[k, j] = i_nodes[j]^k
            end
        end
        Minvt = M1t \ M2t
        for j in eachindex(i_nodes)
            for k in eachindex(i_nodes)
                for n in 1:num_repeats # repeat for the other parameter combinations
                    base_itr_idx = (i + j - 2) * num_repeats + n # exploiting the fact that the pref_itr is last in the product
                    pref_idx = i + k # needs to index the support at i_nodes[k] + lb
                    dvar_idx = (a == pref_alias ? pref_idx : base_itr[base_itr_idx][a] for a in aliases)
                    aug_itr[counter] = (; i = base_itr_idx, c = Minvt[k, j], idx = linear_idxs[dvar_idx...])
                    counter += 1
                end
            end
        end
    end
    dvar = data.infvar_mappings[dref]
    g = p -> p.i => p.c * dvar[p.idx]
    ExaModels.constraint!(core, base_con, Base.Generator(g, aug_itr))
    return 
end
function _add_deriv_equations(core, dref, data, method)
    error("Unsupported derivative reformulation method `$method`.")
end

# Add the approximation equations for each derivative variable
function _add_derivative_approximations(
    core::ExaModels.ExaCore, 
    data::MappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for dref in InfiniteOpt.all_derivatives(inf_model)
        _add_deriv_equations(core, dref, data, InfiniteOpt.derivative_method(dref))
    end
    return
end

# Add the constraints needed for piecewise constant control variables
function _add_collocation_restrictions(
    core::ExaModels.ExaCore, 
    data::MappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for (pidx, vidxs) in inf_model.piecewise_vars
        # gather the basic information
        pref = InfiniteOpt._make_variable_ref(inf_model, pidx)
        if !InfiniteOpt.has_generative_supports(pref)
            continue
        end
        pref_group = InfiniteOpt._object_number(pref)
        pref_alias = data.group_alias[pref_group]
        # make the base pref iterator
        info = InfiniteOpt.generative_support_info(pref)
        num_nodes = length(info.support_basis)
        num_supps = InfiniteOpt.num_supports(pref, label = InfiniteOpt.All)
        ubs = repeat(2+num_nodes:num_nodes+1:num_supps, inner = num_nodes)
        pts = filter(i -> !(i in ubs), 2:num_supps-1)
        pref_itr = [(i1 = ub, i2 = pt) for (ub, pt) in zip(ubs, pts)]
        # make the constraints for each infinite variable
        for vidx in vidxs
            vref = InfiniteOpt._make_variable_ref(inf_model, vidx)
            group_idxs = InfiniteOpt._object_numbers(vref)
            aliases = (data.group_alias[g] for g in group_idxs)
            itrs = (g == pref_group ? pref_itr : data.base_itrs[g] for g in group_idxs)
            itr = vec([merge(i...) for i in Iterators.product(itrs...)])
            ivar = data.infvar_mappings[vref]
            inds1 = Tuple(a == pref_alias ? :i1 : a for a in aliases)
            inds2 = Tuple(a == pref_alias ? :i2 : a for a in aliases)
            f = p -> ivar[values(p[inds1])...] - ivar[values(p[inds2])...]
            ExaModels.constraint(core, Base.Generator(f, itr))
        end
    end
    return
end

# Define warning message for when measure heuristics fail
const _ObjMeasureExpansionWarn = string(
    "Unable to convert objective measures into a form that is ",
    "efficient for ExaModels using existing heuristics. Performance ",
    "may be significantly degraded. Try simplying the objective structure."
)

# Write a finite expression `expr` in a single objective term (this is a generic fallback)
function _add_generic_objective_term(core, expr, data)
    ExaModels.objective(core, _exafy(expr, (;), data), [(;)])
    return
end

# Helper function for generate an iterator based on measure data 
function _make_measure_itr(mdata, data)
    prefs = InfiniteOpt.parameter_refs(mdata) # either a single ref or a vector of refs
    if prefs isa Vector
        supps = eachcol(InfiniteOpt.supports(mdata))
    else
        supps = InfiniteOpt.supports(mdata)
    end
    w = InfiniteOpt.weight_function(mdata)
    coeffs = InfiniteOpt.coefficients(mdata) .* w.(supps)
    prefs = InfiniteOpt.parameter_refs(mdata) # either a single ref or a vector of refs
    group = InfiniteOpt._object_number(first(prefs))
    @assert length(prefs) == length(first(data.base_itrs[group])) - 1 # we don't allow partially measured dependent parameters
    alias = data.group_alias[group]
    aliases = map(p -> data.param_alias[p], prefs)
    return [(; :c => c, alias => data.support_to_index[group, s], zip(aliases, s)...) for (c, s) in zip(coeffs, supps)]
end

# Recursively extract expresion and iterator to be included in the objective
function _process_measure_sum(vref, data, prev_itr = nothing)
    mexpr = InfiniteOpt.measure_function(vref)
    mdata = InfiniteOpt.measure_data(vref)
    curr_itr = _make_measure_itr(mdata, data)
    if isnothing(prev_itr)
        itr = curr_itr
    else
        itr = [(i[1]..., i[2]..., c = i[1].c * i[2].c) for i in Iterators.product(curr_itr, prev_itr)]
    end
    vrefs = InfiniteOpt._all_function_variables(mexpr)
    if all(v.index_type != InfiniteOpt.MeasureIndex for v in vrefs) # single measure without measures inside of it
        return mexpr, itr
    elseif mexpr isa InfiniteOpt.GeneralVariableRef # single nested measure
        return _process_measure_sum(mexpr, data, itr)
    # TODO add more clever heuristics to avoid expanding
    else # fallback for complex nested measures
        inf_model = JuMP.owner_model(vref)
        @warn _ObjMeasureExpansionWarn
        return InfiniteOpt.expand_measures(mexpr, inf_model), itr
    end
end

# Helper function for adding "affine" terms as independent objective terms 
# Note the `coef` doesn't have to be a constant, it can be an expression that doesn't contain measures
function _add_objective_aff_term(core, coef, vref, data)
    return _add_objective_aff_term(core, coef, vref, vref.index_type, data)
end
function _add_objective_aff_term(core, coef, vref, ::Type{InfiniteOpt.MeasureIndex}, data)
    # process the measure structure recursively as needed
    mexpr, itr = _process_measure_sum(vref, data)
    # update the iterator to include parameter function values
    if !isempty(data.pfunc_info)
        itr = _add_parameter_functions_to_itr(mexpr, itr, data)
    end
    # prepare the examodel expression tree
    itr_par = ExaModels.Par(typeof(first(itr)))
    em_expr = itr_par.c * _exafy(coef * mexpr, itr_par, data)
    # add the term to the objective
    ExaModels.objective(core, em_expr, itr)
    return
end
function _add_objective_aff_term(core, coef, vref, _, data)
    _add_generic_objective_term(core, coef * vref, data)
    return
end

# Add the objective from an InfiniteModel to an ExaCore
function _add_objective(
    core::ExaModels.ExaCore,
    expr::JuMP.AbstractJuMPScalar,
    data::MappingData, 
    inf_model::InfiniteOpt.InfiniteModel
    )
    vrefs = InfiniteOpt._all_function_variables(expr)
    if any(v.index_type == InfiniteOpt.MeasureIndex for v in vrefs)
        @warn _ObjMeasureExpansionWarn
    end
    new_expr = InfiniteOpt.expand_measures(expr, inf_model)
    _add_generic_objective_term(core, new_expr, data)
    return
end
function _add_objective(
    core::ExaModels.ExaCore,
    vref::InfiniteOpt.GeneralVariableRef, # can be finite var, point var, finite param, or measure that fully evaluates the measures inside
    data::MappingData, 
    ::InfiniteOpt.InfiniteModel
    )
    _add_objective_aff_term(core, 1.0, vref, data)
    return
end
function _add_objective(
    core::ExaModels.ExaCore,
    aff::JuMP.GenericAffExpr,
    data::MappingData,
    ::InfiniteOpt.InfiniteModel
    )
    # TODO should we check if there are a lot of terms?
    for (coef, vref) in JuMP.linear_terms(aff)
        _add_objective_aff_term(core, coef, vref, data)
    end
    c = JuMP.constant(aff)
    if !iszero(c)
        ExaModels.objective(core, ExaModel.Null(c))
    end
    return
end
function _add_objective(
    core::ExaModels.ExaCore,
    quad::InfiniteOpt.GenericQuadExpr,
    data::MappingData, 
    inf_model::InfiniteOpt.InfiniteModel
    )
    # process the quadratic terms
    for (coef, vref1, vref2) in JuMP.quad_terms(quad)
        if vref1.index_type == InfiniteOpt.MeasureIndex && vref2.index_type == InfiniteOpt.MeasureIndex
            # TODO see if we can avoid the generic fallback in this case
            @warn _ObjMeasureExpansionWarn
            new_expr = InfiniteOpt.expand_measures(coef * vref1 * vref2, inf_model)
            _add_generic_objective_term(core, new_expr, data)
        elseif vref1.index_type == InfiniteOpt.MeasureIndex
            _add_objective_aff_term(core, coef * vref2, vref1, data)
        else
            _add_objective_aff_term(core, coef * vref1, vref2, data)
        end
    end
    # add the affine terms
    _add_objective(core, quad.aff, data, inf_model)
    return
end
# TODO add heuristics for nonlinear expressions

# Fill an ExaCore and MappingData with an InfiniteModel
function build_exa_core!(
    core::ExaModels.ExaCore, 
    data::MappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    # initial setup
    _build_base_iterators(data, inf_model)
    # add the variables and appropriate mappings
    _add_finite_variables(core, data, inf_model)
    _add_infinite_variables(core, data, inf_model) # includes derivatives
    _add_parameter_functions(data, inf_model)
    _add_semi_infinite_variables(data, inf_model)
    _add_point_variables(data, inf_model)
    # add the constraints
    _add_constraints(core, data, inf_model)
    _add_derivative_approximations(core, data, inf_model)
    _add_collocation_restrictions(core, data, inf_model)
    # add the objective if we have one
    expr = JuMP.objective_function(inf_model)
    sense = JuMP.objective_sense(inf_model)
    if sense == _MOI.MAX_SENSE
        _add_objective(core, -expr, data, inf_model)
    elseif sense == _MOI.MIN_SENSE
        _add_objective(core, expr, data, inf_model)
    end
    return
end

# Create an ExaModel from an InfiniteModel
function exa_model(inf_model::InfiniteOpt.InfiniteModel; backend = nothing)
    core = ExaModels.ExaCore(; backend = backend)
    data = MappingData()
    build_exa_core!(core, data, inf_model)
    return ExaModels.ExaModel(core), data 
end

end # end of module