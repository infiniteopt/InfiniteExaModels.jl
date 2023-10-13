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
            Dict{InfiniteOpt.GeneralVariableRef, Tuple{ExaModels.Variable, Vector{Any}}}()
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
                supp = [s for nt in i for s in Iterators.drop(values(nt), 2)]
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
    # return the desired metadata
    return data.infvar_mappings[ivref], indexing
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
        data.semivar_info[vref] = _process_semi_infinite_var(vref, data)
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
# Note: these all will have access to input `p` which is a named tuple from the iterator
function _map_variable(vref, data)
    return _map_variable(vref, vref.index_type, data)
end
function _map_variable(
    vref, 
    ::Type{V}, 
    data
    ) where V <: Union{InfiniteOpt.FiniteVariableIndex}
    return data.finvar_mappings[vref]
end
function _map_variable(
    vref, 
    ::Type{V}, 
    data
    ) where V <: Union{InfiniteOpt.PointVariableIndex}
    if haskey(data.finvar_mappings, vref) 
        return data.finvar_mappings[vref]
    else
        var = _process_point_var(vref, data)
        data.finvar_mappings[vref] = var # TODO consider whether we want to save this or not
        return var
    end
end
function _map_variable(
    vref, 
    ::Type{V}, 
    data
    ) where V <: Union{InfiniteOpt.InfiniteVariableIndex, InfiniteOpt.DerivativeIndex}
    group_idxs = InfiniteOpt._object_numbers(vref)
    ex = Expr(:ref, :($(data.infvar_mappings[vref])))
    append!(ex.args, :(p.$(data.group_alias[i])) for i in group_idxs)
    return ex
end
function _map_variable(
    vref, 
    ::Type{V}, 
    data
    ) where V <: InfiniteOpt.SemiInfiniteVariableIndex
    if haskey(data.semivar_info, vref)
        ivar, inds = data.semivar_info[vref]
    else
        ivar, inds = _process_semi_infinite_var(vref, data)
        data.semivar_info[vref] = (ivar, inds) # TODO consider whether we should save this
    end
    ex = Expr(:ref, :($ivar))
    append!(ex.args, i isa Int ? i : :(p.$i) for i in inds)
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

# Make an a julia function for generators of the form `p -> expr` where `p` is a namedtuple
function _make_expr_ast(expr, data)
    # TODO add operator mapper argument
    return InfiniteOpt.map_expression_to_ast(v -> _map_variable(v, data), expr)
end

# Make generator function from AST expression
function _make_gen_func(expr_code)
    func_code = :(p -> $expr_code) # p is a namedtuple with all the parameters
    return x -> Base.invokelatest(eval(func_code), x)
end

# Construct a Julia function that is compatible with a generator for building ExaModels
function _make_expr_gen_func(expr, data)
    return _make_gen_func(_make_expr_ast(expr, data))
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
        if isempty(inf_model.constraints[JuMP.index(cref)].measure_indices)
            expr = JuMP.jump_function(constr)
        else
            @warn "Constrained measures can lead to poor performance with ExaModels."
            expr = InfiniteOpt.expand_measures(JuMP.jump_function(constr), inf_model)
        end
        set = JuMP.moi_set(constr)
        group_idxs = InfiniteOpt._object_numbers(cref)
        # make the expression generator
        func = _make_expr_gen_func(expr, data)
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
            error("`DomainRestrictions` are not currently supported by ExaModels.")
        end
        gen = Base.Generator(func, itr)
        # get the constraint bounds
        lb, ub = _get_constr_bounds(set)
        # create the ExaModels constraint
        con = ExaModels.constraint(core, gen, lcon = lb, ucon = ub)
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
        inds1 = (i == g_alias ? :i1 : i for i in inds)
        inds2 = (i == g_alias ? :i2 : i for i in inds)
        f = p -> dvar[(p[i] for i in d_inds)...] * p.dt - ivar[(i isa Int ? i : p[i] for i in inds1)...] + ivar[(i isa Int ? i : p[i] for i in inds2)...]
    elseif vref.index_type == InfiniteOpt.InfiniteVariableIndex
        ivar = data.infvar_mappings[vref]
        inds1 = (g == pref_group ? :i1 : data.group_alias[g] for g in group_idxs)
        inds2 = (g == pref_group ? :i2 : data.group_alias[g] for g in group_idxs)
        f = p -> dvar[(p[i] for i in d_inds)...] * p.dt - ivar[(p[i] for i in inds1)...] + ivar[(p[i] for i in inds2)...]
    else # TODO try to make measures work...
        error("Derivatives that act on references with index `$(vref.index_type)` are not supported.")
    end
    # add the equations to core
    gen = Base.Generator(f, itr)
    ExaModels.constraint(core, gen) # TODO add mapping?
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

# Define warning message for when measure heuristics fail
const _ObjMeasureExpansionWarn = string(
    "Unable to convert objective measures into a form that is",
    "efficient for ExaModels using existing heuristics. Performance",
    "may be significantly degraded. Try simplying the objective structure."
)

# Write a finite expression `expr` in a single objective term (this is a generic fallback)
function _add_generic_objective_term(core, expr, data)
    func = _make_expr_gen_func(expr, data)
    ExaModels.objective(core, Base.Generator(func, [(;)]))
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
# TODO finish to function to work properly with recursion for layers of nested measures
function _process_measure_sum(vref, data)
    mexpr = InfiniteOpt.measure_function(vref)
    mdata = InfiniteOpt.measure_data(vref)
    itr = _make_measure_itr(mdata, data)
    vrefs = InfiniteOpt._all_function_variables(mexpr)
    if all(v.index_type != InfiniteOpt.MeasureIndex for v in vrefs) # single measure without measures inside of it
        return mexpr, itr
    elseif mexpr == InfiniteOpt.GeneralVariableRef # single nested measure
        #TODO finish
        error("Nested measures not finished.")
    else # fallback for complex nested measures
        inf_model = JuMP.owner_model(vref)
        @warn _ObjMeasureExpansionWarn
        return InfiniteOpt.expand_measures(mexpr, inf_model), itr
    end
end

# Helper function for adding affine terms as independent objective terms
function _add_objective_aff_term(core, coef, vref, data)
    return _add_objective_aff_term(core, coef, vref, vref.index_type, data)
end
function _add_objective_aff_term(core, coef, vref, ::Type{InfiniteOpt.MeasureIndex}, data)
    mexpr, itr = _process_measure_sum(vref, data)
    expr_code = _make_expr_ast(coef * mexpr, data)
    func = _make_gen_func(:(p.c * $expr_code))
    ExaModels.objective(core, Base.Generator(func, itr))
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
        error("ExaModels does not support objectives with constant terms.")
        # ExaModels.objective(core, c for _ in 1:1) # TODO see if we can make this work
    end
    return
end

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
    _add_semi_infinite_variables(data, inf_model)
    _add_point_variables(data, inf_model)
    # add the constraints
    _add_constraints(core, data, inf_model)
    _add_derivative_approximations(core, data, inf_model)
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
function exa_model(inf_model::InfiniteOpt.InfiniteModel)
    core = ExaModels.ExaCore()
    data = MappingData()
    build_exa_core!(core, data, inf_model)
    return ExaModels.ExaModel(core), data 
end

end # end of module