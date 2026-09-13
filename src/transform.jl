# Create the support iterators for each infinite parameter group and add to the mapping data
function _build_base_iterators(
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    # gather the individual infinite parameter groups
    prefs = InfiniteOpt.parameter_refs(inf_model)
    # build the iterator for each group of infinite parameters
    for group in prefs # group will either be singular ref or a vector of refs
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
        supp_dict = InfiniteOpt.core_object(first(group)).supports
        supps = keys(supp_dict)
        labels = [supp_dict[s] for s in supps]
        group_idx = length(data.group_alias)
        for (i, s) in enumerate(supps)
            data.support_to_index[group_idx, s] = i
        end
        # create the iterator which is a vector of named tuples
        itr = [(; itr_sym => i, zip(aliases, s)...) for (i, s) in enumerate(supps)]
        # add the iterator to `data` and other helpful metadata
        push!(data.base_itrs, itr)
        push!(data.support_labels, labels)
        push!(data.has_internal_supps, InfiniteOpt.has_internal_supports(first(group)))
    end
    return
end

# Ensure the variable is continuous
function _ensure_continuous(info)
    if info.binary || info.integer
        error("Integer variables are not supported by ExaModels.")
    end
end

# Process info value
_process_value(val::Real, supp) = val
_process_value(pf::InfiniteOpt.ParameterFunction, supp) = pf(supp)

## Determine the bounds of an InfiniteOpt variable
# Real bound and start value
function _get_variable_bounds_and_start(
    info::JuMP.VariableInfo{<:Real, <:Real, <:Real, <:Real}
    )
    lb = -Inf
    ub = Inf
    start = 0.0
    if info.has_fix
        lb = info.fixed_value
        ub = lb
    end
    if info.has_lb
        lb = info.lower_bound
    end
    if info.has_ub
        ub = info.upper_bound
    end
    if info.has_start
        start = info.start
    end
    return lb, ub, start
end
# At least one bound/start is a function
function _get_variable_bounds_and_start(info::JuMP.VariableInfo, itrs)
    # set up the collection arrays
    dims = Tuple(length(itr) for itr in itrs)
    lin_idxs = LinearIndices(dims)
    lb = fill(-Inf, length(lin_idxs))
    ub = fill(Inf, length(lin_idxs))
    start = fill(0.0, length(lin_idxs))
    # iterate over all support combinations and fill in the arrays
    for i in Iterators.product(itrs...)
        supp = [s for nt in i for s in Iterators.drop(values(nt), 1)]
        if info.has_fix
            val = _process_value(info.fixed_value, supp)
            lb[lin_idxs[first.(i)...]] = val
            ub[lin_idxs[first.(i)...]] = val
        end
        if info.has_lb
            lb[lin_idxs[first.(i)...]] = _process_value(info.lower_bound, supp)
        end
        if info.has_ub
            ub[lin_idxs[first.(i)...]] = _process_value(info.upper_bound, supp)
        end
        if info.has_start
            start[lin_idxs[first.(i)...]] = _process_value(info.start, supp)
        end
    end
    return lb, ub, start
end

# Get the name from a GeneralVariableRef
function _get_name(vref::InfiniteOpt.GeneralVariableRef, default_name = "var")
    raw_name = JuMP.name(vref)
    sym_name = isempty(raw_name) ? Symbol("$(default_name)$(vref.raw_index)") : Symbol(raw_name)
    return sym_name
end

# Add all the finite variables from an InfiniteModel to a ExaCore
function _add_finite_variables(
    core::ExaModels.ExaCore, 
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    vrefs = JuMP.all_variables(inf_model, InfiniteOpt.FiniteVariable)
    core, ex_vars = ExaModels.add_var(core, length(vrefs), name = Val(:finvar))
    for (i, vref) in enumerate(vrefs)
        info = InfiniteOpt.core_object(vref).info # JuMP.VariableInfo
        _ensure_continuous(info)
        lb, ub, start = _get_variable_bounds_and_start(info)
        ex_var = ex_vars[i]
        data.finvar_mappings[vref] = ex_var
        core.lvar[ex_var.i] = lb
        core.uvar[ex_var.i] = ub
        core.x0[ex_var.i] = start
        data.var_to_grouped_var[vref] = ex_vars
    end
    return core
end

# Add all the finite parameters from an InfiniteModel to a ExaCore
function _add_finite_parameters(
    core::ExaModels.ExaCore, 
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    prefs = JuMP.all_variables(inf_model, InfiniteOpt.FiniteParameter)
    core, ex_pars = ExaModels.add_par(core, length(prefs))
    offset = ex_pars.offset
    for (i, pref) in enumerate(prefs)
        data.param_mappings[pref] = ExaModels.Parameter((1,), 1, offset + i - 1, nothing)
        core.θ[offset + i] = InfiniteOpt.parameter_value(pref)
        data.var_to_grouped_var[pref] = ex_pars
    end
    return core
end

# Add all the infinite variables (and derivatives) from an InfiniteModel to a ExaCore
function _add_infinite_variables(
    core::ExaModels.ExaCore, 
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    # get the raw variables
    ivrefs = JuMP.all_variables(inf_model, InfiniteOpt.InfiniteVariable)
    InfiniteOpt.reformulate_high_order_derivatives!(inf_model)
    drefs = InfiniteOpt.all_derivatives(inf_model)
    vrefs = append!(ivrefs, drefs)
    # sort the variables by parameter groups
    if length(data.base_itrs) > 1
        group_to_vrefs = Dict{Vector{Int}, Vector{InfiniteOpt.GeneralVariableRef}}()
        for vref in vrefs
            group_idxs = InfiniteOpt.parameter_group_int_indices(vref)
            if !haskey(group_to_vrefs, group_idxs)
                group_to_vrefs[group_idxs] = [vref]
            else
                push!(group_to_vrefs[group_idxs], vref)
            end
        end
    else
        group_idxs = InfiniteOpt.parameter_group_int_indices(first(vrefs))
        group_to_vrefs = Dict(group_idxs => vrefs)
    end
    # add each group of variables in group_to_vrefs to the ExaCore
    for (group_idxs, vrefs) in group_to_vrefs
        itrs = map(i -> data.base_itrs[i], group_idxs)
        dims = Tuple(length(itr) for itr in itrs)
        core, ex_vars = ExaModels.add_var(core, dims..., length(vrefs), name = Val(:infvar))
        offset = ex_vars.offset
        for vref in vrefs
            info = InfiniteOpt.core_object(vref).info # JuMP.VariableInfo
            _ensure_continuous(info)
            lb, ub, start = _get_variable_bounds_and_start(info, itrs)
            vname = _get_name(vref, vref in drefs ? "deriv" : "infvar")
            data.infvar_mappings[vref] = ExaModels.Variable(dims, length(lb), offset, vname, nothing)
            copyto!(@view(core.lvar[offset+1:offset+length(lb)]), lb)
            copyto!(@view(core.uvar[offset+1:offset+length(ub)]), ub)
            copyto!(@view(core.x0[offset+1:offset+length(start)]), start)
            offset += length(lb)
            data.var_to_grouped_var[vref] = ex_vars
        end
    end
    return core
end

# Process all the parameter function from an InfiniteModel and add to an ExaCore
function _add_parameter_functions(
    core::ExaModels.ExaCore,
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )  
    pfrefs = InfiniteOpt.all_parameter_functions(inf_model)
    iszero(length(pfrefs)) && return core 
    # sort the parameter functions by parameter groups
    if length(data.base_itrs) > 1
        group_to_pfrefs = Dict{Vector{Int}, Vector{InfiniteOpt.GeneralVariableRef}}()
        for pfref in pfrefs
            group_idxs = InfiniteOpt.parameter_group_int_indices(pfref)
            if !haskey(group_to_pfrefs, group_idxs)
                group_to_pfrefs[group_idxs] = [pfref]
            else
                push!(group_to_pfrefs[group_idxs], pfref)
            end
        end
    else
        group_idxs = InfiniteOpt.parameter_group_int_indices(first(pfrefs))
        group_to_pfrefs = Dict(group_idxs => pfrefs)
    end
    # add each group of parameter functions to the ExaCore
    for (group_idxs, group_pfrefs) in group_to_pfrefs
        itrs = map(i -> data.base_itrs[i], group_idxs)
        dims = Tuple(length(itr) for itr in itrs)
        core, ex_pars = ExaModels.add_par(core, dims..., length(group_pfrefs))
        offset = ex_pars.offset
        for pfref in group_pfrefs
            pfunc = InfiniteOpt.core_object(pfref)
            lin_idxs = LinearIndices(dims)
            vals = Vector{Float64}(undef, length(lin_idxs))
            for i in Iterators.product(itrs...)
                supp = [s for nt in i for s in Iterators.drop(values(nt), 1)]
                vals[lin_idxs[first.(i)...]] = pfunc(supp)
            end
            copyto!(@view(core.θ[offset+1:offset+length(vals)]), vals)
            data.param_mappings[pfref] = ExaModels.Parameter(dims, length(vals), offset, nothing)
            offset += length(vals)
            data.var_to_grouped_var[pfref] = ex_pars
        end
    end
    return core
end

# Helper function for processing semi-infinite variables
function _process_semi_infinite_var(vref, data)
    # get basic info from InfiniteOpt
    ivref = InfiniteOpt.infinite_variable_ref(vref)
    orig_group_idxs = InfiniteOpt.parameter_group_int_indices(ivref)
    raw_prefs = InfiniteOpt.raw_parameter_refs(ivref) # type `InfiniteOpt.VectorTuple`
    group_idxs = InfiniteOpt.parameter_group_int_indices(vref)
    eval_supp = InfiniteOpt.eval_support(vref)
    # create metadata vector `indexing` 
    indexing = Vector{Any}(undef, length(orig_group_idxs))
    for (i, g) in enumerate(orig_group_idxs)
        if g in group_idxs
            indexing[i] = data.group_alias[g] # get the group alias for indexing
        elseif iszero(raw_prefs.dimensions[i])
            supp = eval_supp[first(raw_prefs.ranges[i])]
            indexing[i] = data.support_to_index[g, supp] # store the support index
        else
            supp = eval_supp[raw_prefs.ranges[i]]
            indexing[i] = data.support_to_index[g, supp] # store the support index
        end
    end
    # store the desired information
    if ivref.index_type == InfiniteOpt.ParameterFunctionIndex
        mapped_var = data.param_mappings[ivref]
    else
        mapped_var = data.infvar_mappings[ivref]
    end
    if haskey(data.var_to_grouped_var, ivref)
        data.var_to_grouped_var[vref] = data.var_to_grouped_var[ivref]
    end
    return data.semivar_info[vref] = (mapped_var, indexing)
end

# Update the bounds and start value of a ExaModels.Var based on RestrictedDomainInfo
function _update_bounds_and_start(core, info, var)
    if info.active_lower_bound_info
        core.lvar[var.i] = isnan(info.lower_bound) ? -Inf : info.lower_bound
    end
    if info.active_upper_bound_info
        core.uvar[var.i] = isnan(info.upper_bound) ? Inf : info.upper_bound
    end
    if info.active_fix_info
        core.lvar[var.i] = isnan(info.fixed_value) ? -Inf : info.fixed_value
        core.uvar[var.i] = isnan(info.fixed_value) ? Inf : info.fixed_value
    end
    if info.active_start_info
        core.x0[var.i] = info.start_value
    end
    return
end

# Add all the semi-infinite variables from an InfiniteModel to a ExaCore
# In other words, create helpful metadata to be used by `_map_variable`
function _add_semi_infinite_variables(
    core::ExaModels.ExaCore,
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for vref in JuMP.all_variables(inf_model, InfiniteOpt.SemiInfiniteVariable)
        # collect the basic information fields from InfiniteOpt and save
        mapped_var, indexing = _process_semi_infinite_var(vref, data)
        # updated the bounds and start value if needed
        info = InfiniteOpt.core_object(vref).info # InfiniteOpt.RestrictedDomainInfo
        if info.active_lower_bound_info || info.active_upper_bound_info ||
           info.active_fix_info || info.active_start_info
            semivar_idxs = (idx isa Int ? idx : 1:mapped_var.size[i] for (i, idx) in enumerate(indexing))   
            semivar_itr = Iterators.product(semivar_idxs...)
            for idx in semivar_itr
                var = mapped_var[idx...]
                _update_bounds_and_start(core, info, var)
            end
        end
    end
    return 
end

# Helper function for processing point variables
function _process_point_var(vref, data)
    ivref = InfiniteOpt.infinite_variable_ref(vref)
    raw_supp = InfiniteOpt.raw_parameter_values(vref)
    prefs = InfiniteOpt.raw_parameter_refs(ivref)
    supp = Tuple(raw_supp, prefs)
    if any(d >= 2 for d in prefs.dimensions)
        supp = Tuple(InfiniteOpt.Collections.vectorize(s)[1] for s in supp)
    end
    group_idxs = InfiniteOpt.parameter_group_int_indices(ivref)
    idxs = Tuple(data.support_to_index[i, s] for (i, s) in zip(group_idxs, supp))
    pt = data.infvar_mappings[ivref][idxs...]
    if haskey(data.var_to_grouped_var, ivref)
        data.var_to_grouped_var[vref] = data.var_to_grouped_var[ivref]
    end
    data.point_indicies[vref] = idxs
    return data.finvar_mappings[vref] = pt
end

# Add all the point variables from an InfiniteModel to a ExaCore
function _add_point_variables(
    core::ExaModels.ExaCore,
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for vref in JuMP.all_variables(inf_model, InfiniteOpt.PointVariable)
        # store the index mapping for the point variable
        pt = _process_point_var(vref, data)
        # update the bounds and start value if needed
        info = InfiniteOpt.core_object(vref).info # InfiniteOpt.RestrictedDomainInfo
        _update_bounds_and_start(core, info, pt)
    end
    return 
end

# Get the index parameters for a variable reference (used by `_map_variable`)
function _index_params(
    vref::InfiniteOpt.GeneralVariableRef,
    data::ExaMappingData
    )
    _index_params(vref, vref.index_type, data)
end
function _index_params(
    vref::InfiniteOpt.GeneralVariableRef,
    ::Type{V},
    data::ExaMappingData
    ) where V <: Union{InfiniteOpt.InfiniteVariableIndex, InfiniteOpt.DerivativeIndex, InfiniteOpt.ParameterFunctionIndex}
    group_idxs = InfiniteOpt.parameter_group_int_indices(vref)
    data_src = ExaModels.DataSource()
    return (data_src[data.group_alias[i]] for i in group_idxs)
end
function _index_params(
    vref::InfiniteOpt.GeneralVariableRef,
    ::Type{InfiniteOpt.SemiInfiniteVariableIndex},
    data::ExaMappingData
    )
    if !haskey(data.semivar_info, vref)
        _process_semi_infinite_var(vref, data)
    end
    _, inds = data.semivar_info[vref]
    data_src = ExaModels.DataSource()
    return (i isa Int ? i : data_src[i] for i in inds)
end
function _index_params(
    vref::InfiniteOpt.GeneralVariableRef,
    ::Type{InfiniteOpt.PointVariableIndex},
    data::ExaMappingData
    )
    if !haskey(data.finvar_mappings, vref)
        _process_point_var(vref, data)
    end
    return data.point_indicies[vref]
end
function _index_params(
    vref::InfiniteOpt.GeneralVariableRef,
    type,
    data::ExaMappingData
    )
    return ()
end 
 
# Add user-defined operators to ExaModels
function _add_user_operators(inf_model::InfiniteOpt.InfiniteModel)
    for op in InfiniteOpt.added_nonlinear_operators(inf_model)
        if haskey(_op_mappings, op.name)
            continue
        elseif op.dim > 1
            error("InfiniteExaModels.jl does not currently support multivariate user-defined nonlinear operators. ",
                  "If you need support for this operator, please open an issue.")
        elseif op.∇²f === nothing || op.∇f === nothing
            error("InfiniteExaModels.jl does not support user-defined nonlinear operators without a gradient and Hessian.")
        else
            @eval begin 
                ExaModels.@register_univariate($(op).f, $(op).∇f, $(op).∇²f) 
            end
            _op_mappings[op.name] = op.f
        end
    end
    return
end

# Map variable references based on their underlying type (used by `_exafy`)
function _map_variable(vref, ::Type{InfiniteOpt.FiniteVariableIndex}, data)
    return data.finvar_mappings[vref]
end
function _map_variable(vref, ::Type{InfiniteOpt.PointVariableIndex}, data)
    if haskey(data.finvar_mappings, vref) 
        return data.finvar_mappings[vref]
    else
        return _process_point_var(vref, data)
    end
end
function _map_variable(
    vref,
    ::Type{V},
    data
    ) where V <: Union{InfiniteOpt.InfiniteVariableIndex, InfiniteOpt.DerivativeIndex}
    idx_pars = _index_params(vref, V, data)
    return data.infvar_mappings[vref][idx_pars...]
end
function _map_variable(vref, ::Type{InfiniteOpt.SemiInfiniteVariableIndex}, data)
    idx_pars = _index_params(vref, InfiniteOpt.SemiInfiniteVariableIndex, data)
    ivar, _ = data.semivar_info[vref]
    return ivar[idx_pars...]
end
function _map_variable(vref, ::Type{<:InfiniteOpt.InfiniteParameterIndex}, data)
    return ExaModels.DataSource()[data.param_alias[vref]]
end
function _map_variable(vref, ::Type{InfiniteOpt.FiniteParameterIndex}, data)
    return data.param_mappings[vref][1]
end
function _map_variable(vref, ::Type{InfiniteOpt.ParameterFunctionIndex}, data)
    idx_pars = _index_params(vref, InfiniteOpt.ParameterFunctionIndex, data)
    return data.param_mappings[vref][idx_pars...]
end
function _map_variable(vref, IdxType, data)
    error("Unable to add `$vref` to an ExaModel, it's index type `$IdxType`" *
          " is not yet supported by InfiniteExaModels.")
end

# Convert as InfiniteOpt expression into a ExaModel expression
function _exafy(vref::InfiniteOpt.GeneralVariableRef, data)
    return _map_variable(vref, vref.index_type, data)
end
function _exafy(c::Number, data)
    return c
end
function _exafy(
    aff::JuMP.GenericAffExpr{C, InfiniteOpt.GeneralVariableRef},
    data
    ) where {C}
    c = JuMP.constant(aff)
    if !isempty(aff.terms)
        ex = sum(begin
            v_ex = _exafy(v, data)
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
    data
    ) where {C}
    aff = _exafy(quad.aff, data)
    if !isempty(quad.terms)
        ex = sum(begin 
            if v1 == v2
                v_ex = _exafy(v1, data) 
                isone(c) ? abs2(v_ex) : c * abs2(v_ex)
            else
                v1_ex = _exafy(v1, data) 
                v2_ex = _exafy(v2, data) 
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
    data
    )
    return _nl_op(nl.head)((_exafy(a, data) for a in nl.args)...)
end

# Check if NamedTuple iterator respects the restriction
function _support_in_restriction(restriction, itr, data)
    supp = [itr[data.param_alias[p]] for p in restriction.parameter_refs]
    return restriction(supp)
end

# Prepare the constraint iterator over the infinite parameters
function _get_constraint_iterator(cref, data)
    group_idxs = InfiniteOpt.parameter_group_int_indices(cref)
    # prepare the iterator of NamedTuples (contains support values, iterator values, and constants from parameter functions)
    if isempty(group_idxs) # we have a finite constraint
        itr = [(;)]
    elseif length(group_idxs) == 1 # we only depend on one independent infinite parameter
        itr = data.base_itrs[first(group_idxs)]
    else # we depend on multiple independent infinite parameters
        itrs = map(i -> data.base_itrs[i], group_idxs)
        itr = vec([merge(i...) for i in Iterators.product(itrs...)])
    end
    # Remove any elements of the iterator that violate the domain restriction
    if InfiniteOpt.has_domain_restriction(cref)
        restriction = InfiniteOpt.domain_restriction(cref)
        itr = filter(i -> _support_in_restriction(restriction, i, data), itr)
    end
    return itr
end

# Finalize exafied expressions to avoid scalars
_finalize_expr(expr) = expr
_finalize_expr(c::Real) = ExaModels.Null(c)

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
    error("Constraint set `$set` is not supported by InfiniteExaModels, ",
          "if you need support for this constraint type, please open an issue.")
end

# Add all the constraints from an InfiniteModel to an ExaCore
function _add_constraints(
    core::ExaModels.ExaCore, 
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for cref in JuMP.all_constraints(inf_model)
        # skip if the constraint is a variable bound or already added (as a grouped constraint)
        InfiniteOpt.is_variable_domain_constraint(cref) && continue
        haskey(data.constraint_mappings, cref) && continue
        # parse the basic information
        constr = JuMP.constraint_object(cref)
        if isempty(inf_model.constraints[JuMP.index(cref)].measure_indices)
            expr = JuMP.jump_function(constr)
        else
            @warn "Constrained measures can lead to poor performance with ExaModels."
            expr = InfiniteOpt.expand_measures(JuMP.jump_function(constr), inf_model)
        end
        set = JuMP.moi_set(constr)
       # prepare the constraint iterator
        itr = _get_constraint_iterator(cref, data)
        # create the ExaModels expression tree based on expr
        em_expr = _finalize_expr(_exafy(expr, data))
        # get the constraint bounds
        lb, ub = _get_constr_bounds(set)
        # create the ExaModels constraint
        core, con = ExaModels.add_con(core, em_expr, itr, lcon = lb, ucon = ub)
        data.constraint_mappings[cref] = con
    end
    return core
end

# Make dispatch type to pass the data needed by `make_reduced_expr`
struct _DerivReductionBackendInfo <: InfiniteOpt.AbstractTransformationBackend
    data::ExaMappingData
    alias_map::Union{Nothing, Dict{InfiniteOpt.GeneralVariableRef, Symbol}}
end

# Extend make_reduced_expr to create an ExaModel expression
function InfiniteOpt.make_reduced_expr(
    vref::InfiniteOpt.GeneralVariableRef,
    pref::InfiniteOpt.GeneralVariableRef,
    supps::Vector{Float64},
    idx,
    dispatch_data::_DerivReductionBackendInfo
    )
    group_idx = InfiniteOpt.parameter_group_int_index(pref)
    data = dispatch_data.data
    data_src = ExaModels.DataSource()
    alias = data.group_alias[group_idx]
    if vref.index_type == InfiniteOpt.SemiInfiniteVariableIndex
        @assert haskey(data.semivar_info, vref)
        ivar, inds = data.semivar_info[vref]
        idx_pars = (begin 
            if i isa Int
                i 
            elseif i == alias
                idx
            else
                data_src[i] 
            end
            end for i in inds)
        grouped_var = data.var_to_grouped_var[vref]
        return grouped_var[idx_pars..., data_src[dispatch_data.alias_map[vref]]]
    else # either an infinite variable or a derivative variable
        group_idxs = InfiniteOpt.parameter_group_int_indices(vref)
        idx_pars = (begin 
            g_alias = data.group_alias[i]
            if g_alias == alias
                idx
            else
                data_src[g_alias]
            end 
            end for i in group_idxs)
        grouped_var = data.var_to_grouped_var[vref]
        return grouped_var[idx_pars..., data_src[dispatch_data.alias_map[vref]]]
    end
end

# Add the approximation equations for each derivative variable
function _add_derivative_approximations(
    core::ExaModels.ExaCore, 
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    # group all the derivatives of the same order, method, and infinite parameter dependencies
    signature_to_derivs = Dict{
        Tuple{InfiniteOpt.GeneralVariableRef, Int, Vector{Int}, DataType}, 
        Tuple{Vector{InfiniteOpt.GeneralVariableRef}, Vector{InfiniteOpt.GeneralVariableRef}}
    }()
    for dref in InfiniteOpt.all_derivatives(inf_model)
        vref = InfiniteOpt.derivative_argument(dref) 
        pref = InfiniteOpt.operator_parameter(dref)
        order = InfiniteOpt.derivative_order(dref)
        group_idxs = InfiniteOpt.parameter_group_int_indices(vref)
        if !haskey(signature_to_derivs, (pref, order, group_idxs, vref.index_type))
            signature_to_derivs[pref, order, group_idxs, vref.index_type] = 
            (InfiniteOpt.GeneralVariableRef[], InfiniteOpt.GeneralVariableRef[])
        end
        push!(signature_to_derivs[pref, order, group_idxs, vref.index_type][1], dref)
        push!(signature_to_derivs[pref, order, group_idxs, vref.index_type][2], vref)
    end
    # iterate over each group of derivatives and add the approximation equations
    for ((pref, order, group_idxs, _), (drefs, vrefs)) in signature_to_derivs
        # gather basic info
        method = InfiniteOpt.derivative_method(drefs[1])
        pref_group = InfiniteOpt.parameter_group_int_index(pref)
        # sort the base support iterator 
        p_alias = data.param_alias[pref]
        base_itr = data.base_itrs[pref_group]
        if length(first(base_itr)) == 2 # we have a independent parameter (the supports are sorted already)
            srt_itr = base_itr
        else # we have a dependent parameter (supports are not sorted)
            srt_itr = sort(base_itr, by = p -> p[p_alias])
        end
        # collect the expression data
        supps = map(p -> p[p_alias], srt_itr)
        idxs, arg_itrs... = InfiniteOpt.derivative_expr_data(drefs[1], order, supps, method)
        # make the iterator
        aliases = Tuple(Symbol("d_arg$i") for i in eachindex(arg_itrs))
        pref_itr = [(; srt_itr[i]..., zip(aliases, args)...) for (i, args...) in zip(idxs, arg_itrs...)]
        itrs = Any[g == pref_group ? pref_itr : data.base_itrs[g] for g in group_idxs]
        push!(
            itrs, 
            [(; :grouped_didx => _get_grouped_idx(dref, data), 
                :grouped_vidx => _get_grouped_idx(vref, data)) 
                for (dref, vref) in zip(drefs, vrefs)
            ])
        itr = length(itrs) > 1 ? vec([merge(i...) for i in Iterators.product(itrs...)]) : pref_itr
        # make the ExaModel expression tree and add the constraint(s)
        data_src = ExaModels.DataSource()
        alias_map = Dict(drefs[1] => :grouped_didx, vrefs[1] => :grouped_vidx)
        em_expr = InfiniteOpt.make_indexed_derivative_expr(
            drefs[1], 
            vrefs[1],
            pref,
            order,
            data_src[data.group_alias[pref_group]],
            supps,
            _DerivReductionBackendInfo(data, alias_map),
            method,
            (data_src[a] for a in aliases)...
        )
        core, _ = ExaModels.add_con(core, em_expr, itr)
    end
    return core
end

# Add the constraints needed for piecewise constant control variables
function _add_collocation_restrictions(
    core::ExaModels.ExaCore, 
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for (pidx, vidxs) in inf_model.piecewise_vars
        # gather the basic information
        pref = InfiniteOpt.GeneralVariableRef(inf_model, pidx)
        if !InfiniteOpt.has_generative_supports(pref)
            continue
        end
        pref_group = InfiniteOpt.parameter_group_int_index(pref)
        pref_alias = data.group_alias[pref_group]
        # make the base pref iterator
        info = InfiniteOpt.generative_support_info(pref)
        num_nodes = length(info.support_basis)
        num_supps = InfiniteOpt.num_supports(pref, label = InfiniteOpt.All)
        ubs = repeat(2+num_nodes:num_nodes+1:num_supps, inner = num_nodes)
        pts = filter(i -> !(i in ubs), 2:num_supps-1)
        pref_itr = [(i1 = ub, i2 = pt) for (ub, pt) in zip(ubs, pts)]
        # group the variables by their input infinite parameters
        group_idxs_to_vrefs = Dict{Vector{Int}, Vector{InfiniteOpt.GeneralVariableRef}}()
        for vidx in vidxs
            vref = InfiniteOpt.GeneralVariableRef(inf_model, vidx)
            group_idxs = InfiniteOpt.parameter_group_int_indices(vref)
            if !haskey(group_idxs_to_vrefs, group_idxs)
                group_idxs_to_vrefs[group_idxs] = InfiniteOpt.GeneralVariableRef[]
            end
            push!(group_idxs_to_vrefs[group_idxs], vref)
        end
        # add the constraints for each group of variables
        for (group_idxs, vrefs) in group_idxs_to_vrefs
            # prepare the iterator
            aliases = (data.group_alias[g] for g in group_idxs)
            itrs = (g == pref_group ? pref_itr : data.base_itrs[g] for g in group_idxs)
            finite_itr = [(; :grouped_vidx => _get_grouped_idx(vref, data)) for vref in vrefs]
            itr = vec([merge(i...) for i in Iterators.product(itrs..., finite_itr)])
            # prepare the variable indices
            data_src = ExaModels.DataSource()
            alias_tuple = (aliases..., :grouped_vidx)
            idx_pars1 = (a == pref_alias ? data_src[:i1] : data_src[a] for a in alias_tuple)
            idx_pars2 = (a == pref_alias ? data_src[:i2] : data_src[a] for a in alias_tuple)
            # create the ExaModel expression tree and add the constraint
            grouped_var = data.var_to_grouped_var[vrefs[1]]
            em_expr = grouped_var[idx_pars1...] - grouped_var[idx_pars2...]
            core, _ = ExaModels.add_con(core, em_expr, itr)
        end
    end
    return core
end

# Define warning message for when measure heuristics fail
const _ObjMeasureExpansionWarn = string(
    "Unable to convert objective measures into a form that is ",
    "efficient for ExaModels using existing heuristics. Performance ",
    "may be significantly degraded. Try simplying the objective structure. ",
    "if you think this form should be supported, please open an issue."
)

# Write a finite expression `expr` in a single objective term (this is a generic fallback)
function _add_generic_objective_term(core, expr, data)
    em_expr = _finalize_expr(_exafy(expr, data))
    return ExaModels.add_obj(core, em_expr, [(;)])[1]
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
    group = InfiniteOpt.parameter_group_int_index(first(prefs))
    @assert length(prefs) == length(first(data.base_itrs[group])) - 1 # we don't allow partially measured dependent parameters
    alias = data.group_alias[group]
    aliases = map(p -> data.param_alias[p], prefs)
    return [(; :c => c, alias => data.support_to_index[group, s], zip(aliases, s)...) for (c, s) in zip(coeffs, supps)]
end

# Does an expression contain vref? (TODO: make this more efficient)
function _has_variable(expr, vref)
    return vref in InfiniteOpt.all_expression_variables(expr)
end

# Determine if expression can be reformulated such that all terms move inside the inner measure
# Assumes that `expr` only contains 1 unique measure
_terms_can_be_moved_inside_measure(::InfiniteOpt.GeneralVariableRef, mref) = true
_terms_can_be_moved_inside_measure(::JuMP.GenericAffExpr, mref) = true
function _terms_can_be_moved_inside_measure(quad::JuMP.GenericQuadExpr, mref)
    return !haskey(quad.terms, JuMP.UnorderedPair(mref, mref))
end
function _terms_can_be_moved_inside_measure(nlp::JuMP.GenericNonlinearExpr, mref)
    m_inds = findall(ex -> _has_variable(ex, mref), nlp.args)
    if nlp.head in (:+, :-)
        return all(_terms_can_be_moved_inside_measure(ex, mref) for ex in nlp.args[m_inds])
    elseif nlp.head == :*
        if length(m_inds) > 1
            return false
        else
            return _terms_can_be_moved_inside_measure(nlp.args[only(m_inds)], mref)
        end
    else
        return false
    end
end

# Recursively extract expression and iterator to be included in the objective
function _process_measure_sum(vref, data, prev_itr = nothing)
    mexpr = InfiniteOpt.measure_function(vref)
    mdata = InfiniteOpt.measure_data(vref)
    curr_itr = _make_measure_itr(mdata, data)
    if isnothing(prev_itr)
        itr = curr_itr
    else
        itr = [(i[1]..., i[2]..., c = i[1].c * i[2].c) for i in Iterators.product(curr_itr, prev_itr)]
    end
    vrefs = InfiniteOpt.all_expression_variables(mexpr)
    mrefs = filter(v -> v.index_type == InfiniteOpt.MeasureIndex, vrefs)
    if isempty(mrefs) # single measure without measures inside of it
        return mexpr, itr
    elseif isone(length(mrefs)) && _terms_can_be_moved_inside_measure(mexpr, only(mrefs))
        mref = only(mrefs)
        inner_mexpr, new_itr = _process_measure_sum(mref, data, itr)
        return InfiniteOpt.map_expression(v -> isequal(v, mref) ? inner_mexpr : v, mexpr), new_itr
    # TODO add more clever heuristics to avoid expanding
    else # fallback for complex nested measures
        inf_model = JuMP.owner_model(vref)
        @warn _ObjMeasureExpansionWarn
        return InfiniteOpt.expand_measures(mexpr, inf_model), itr
    end
end

# Helper function for adding "affine" terms as independent objective terms 
# Note the `coef` doesn't have to be a constant, it can be an expression that doesn't contain measures
function _add_objective_aff_term(core, coef, vref, data, group_repeated_sums = false)
    return _add_objective_aff_term(core, coef, vref, vref.index_type, data, group_repeated_sums)
end
function _add_objective_aff_term(core, coef, vref, ::Type{InfiniteOpt.MeasureIndex}, data, group_repeated_sums)
    # process the measure structure recursively as needed
    mexpr, itr = _process_measure_sum(vref, data)
    # form the exafied expression and iterator
    c = ExaModels.DataSource()[:c]
    if group_repeated_sums
        exafied_expr, finite_itr = _process_candidate_sum_group(mexpr, data)
        if length(finite_itr) > 1
            @info "Successfully grouped $(length(finite_itr)) finite terms together into a single objective pattern."
            final_itr = vec([merge(i...) for i in Iterators.product(itr, finite_itr)])
        else
            final_itr = itr
        end
    else
        exafied_expr = _exafy(mexpr, data)
        final_itr = itr
    end
    # prepare the examodel expression tree
    em_expr = isone(coef) ? c * exafied_expr : _exafy(coef, data) * (c * exafied_expr)
    # add the term to the objective
    core, _ = ExaModels.add_obj(core, _finalize_expr(em_expr), final_itr)
    return core
end
function _add_objective_aff_term(core, coef, vref, _, data, group_repeated_sums)
    expr = isone(coef) ? vref : coef * vref
    return _add_generic_objective_term(core, expr, data)
end

# Add the objective from an InfiniteModel to an ExaCore
function _add_objective(
    core::ExaModels.ExaCore,
    expr::JuMP.AbstractJuMPScalar, # generic fallback (heuristics fail to find a summed measure structure)
    data::ExaMappingData, 
    inf_model::InfiniteOpt.InfiniteModel;
    group_repeated_sums::Bool = false
    )
    vrefs = InfiniteOpt.all_expression_variables(expr)
    if any(v.index_type == InfiniteOpt.MeasureIndex for v in vrefs)
        @warn _ObjMeasureExpansionWarn
    end
    new_expr = InfiniteOpt.expand_measures(expr, inf_model)
    return _add_generic_objective_term(core, new_expr, data)
end
function _add_objective(
    core::ExaModels.ExaCore,
    vref::InfiniteOpt.GeneralVariableRef, # can be finite var, point var, finite param, or measure that fully evaluates the measures inside
    data::ExaMappingData, 
    ::InfiniteOpt.InfiniteModel;
    group_repeated_sums::Bool = false
    )
    return _add_objective_aff_term(core, 1.0, vref, data, group_repeated_sums)
end
function _add_objective(
    core::ExaModels.ExaCore,
    aff::JuMP.GenericAffExpr,
    data::ExaMappingData,
    ::InfiniteOpt.InfiniteModel;
    group_repeated_sums::Bool = false
    )
    # TODO should we check if there are a lot of terms? (use group_repeated_sums)
    for (coef, vref) in JuMP.linear_terms(aff)
        core = _add_objective_aff_term(core, coef, vref, data, group_repeated_sums)
    end
    c = JuMP.constant(aff)
    if !iszero(c)
        core, _ = ExaModels.add_obj(core, ExaModels.Null(c))
    end
    return core
end
function _add_objective(
    core::ExaModels.ExaCore,
    quad::InfiniteOpt.GenericQuadExpr,
    data::ExaMappingData, 
    inf_model::InfiniteOpt.InfiniteModel;
    group_repeated_sums::Bool = false
    )
    # process the quadratic terms
    for (coef, vref1, vref2) in JuMP.quad_terms(quad)
        if vref1.index_type == InfiniteOpt.MeasureIndex && vref2.index_type == InfiniteOpt.MeasureIndex
            # TODO see if we can avoid the generic fallback in this case
            @warn _ObjMeasureExpansionWarn
            new_expr = InfiniteOpt.expand_measures(coef * vref1 * vref2, inf_model)
            core = _add_generic_objective_term(core, new_expr, data)
        elseif vref1.index_type == InfiniteOpt.MeasureIndex
            core = _add_objective_aff_term(core, coef * vref2, vref1, data, group_repeated_sums)
        else
            core = _add_objective_aff_term(core, coef * vref1, vref2, data, group_repeated_sums)
        end
    end
    # add the affine terms
    core = _add_objective(core, quad.aff, data, inf_model)
    return core
end
# TODO add heuristics for nonlinear expressions

# Fill an ExaCore and ExaMappingData with an InfiniteModel
function build_exa_core!(
    core::ExaModels.ExaCore, 
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel;
    group_repeated_algebraic_patterns = false
    )
    # initial setup
    _build_base_iterators(data, inf_model)
    # add the variables and appropriate mappings
    core = _add_finite_parameters(core, data, inf_model)
    core = _add_finite_variables(core, data, inf_model)
    core = _add_infinite_variables(core, data, inf_model)
    core = _add_parameter_functions(core, data, inf_model)
    _add_semi_infinite_variables(core, data, inf_model)
    _add_point_variables(core, data, inf_model)
    # account for user-defined nonlinear operators
    _add_user_operators(inf_model)
    # add the constraints
    if group_repeated_algebraic_patterns
        core = _group_and_add_constraints(core, data, inf_model) # TODO: can eventually replace `_add_constraints` if it works well
        num_grouped_constraints = length(core.cons)
    end
    core = _add_constraints(core, data, inf_model)
    if group_repeated_algebraic_patterns
        num_ungrouped_constraints = length(core.cons) - num_grouped_constraints
    end
    core = _add_derivative_approximations(core, data, inf_model)
    core = _add_collocation_restrictions(core, data, inf_model)
    # add the objective if there is one
    expr = JuMP.objective_function(inf_model)
    sense = JuMP.objective_sense(inf_model)
    if sense != _MOI.FEASIBILITY_SENSE
        core = _add_objective(
            core, 
            expr, 
            data, 
            inf_model, 
            group_repeated_sums = group_repeated_algebraic_patterns
        )
    end
    if group_repeated_algebraic_patterns
        num_con_patterns = length(core.cons)
        num_grouped_constraints = num_con_patterns - num_ungrouped_constraints
        @info "In total, $num_con_patterns constraint pattern(s) was/were added of which $num_grouped_constraints are grouped constraints."
        @info "In total, $(length(core.obj)) objective sum pattern(s) was/were added. Check the logs to determine how many were grouped."
    end
    return core
end

"""
    ExaModels.ExaCore(
        inf_model::InfiniteOpt.InfiniteModel,
        data::ExaMappingData;
        [backend = nothing,
        concrete_core::Bool = false,
        group_repeated_algebraic_patterns = false] # experimental
    )::ExaModels.ExaCore

Create `ExaModels.ExaCore` from `inf_model` using the provided
`ExaMappingData` to store the variable and constraint mappings. 
The setting `concrete_core = true` will create a concrete 
`ExaModels.ExaCore` type, which is useful for performance in some cases.
Optionally, try to aggregate common algebraic constraint and objective patterns
by setting `group_repeated_algebraic_patterns = true`. This is an 
experimental feature that may encounter issues and may be removed/modified in the future.
"""
function ExaModels.ExaCore(
    inf_model::InfiniteOpt.InfiniteModel,
    data::ExaMappingData;
    backend = nothing,
    concrete_core::Bool = false,
    group_repeated_algebraic_patterns = false
    )
    # TODO add support for other float types once InfiniteOpt does
    minimize = JuMP.objective_sense(inf_model) == _MOI.MIN_SENSE
    core = ExaModels.ExaCore(; backend = backend, minimize = minimize, concrete = Val(concrete_core))
    return build_exa_core!(
        core,
        data, 
        inf_model; group_repeated_algebraic_patterns = group_repeated_algebraic_patterns
    )
end

"""
    ExaModels.ExaModel(
        inf_model::InfiniteOpt.InfiniteModel,
        [data::ExaMappingData];
        [backend = nothing,
        concrete_core::Bool = false,
        group_repeated_algebraic_patterns = false] # experimental
    )::ExaModels.ExaModel

Create an `ExaModels.ExaModel` from `inf_model` and store the mappings in
`data`. If `data` is not provided, the mappings cannot be readily extracted.
The `concrete_core` setting will create a concrete `ExaModels.ExaCore` type, 
which is useful for performance in some cases.
Optionally, try to aggregate common algebraic constraint/objective patterns
by setting `group_repeated_algebraic_patterns = true`. This is an 
experimental feature that may encounter issues and may be removed/modified in the future.
"""
function ExaModels.ExaModel(
    inf_model::InfiniteOpt.InfiniteModel,
    data::ExaMappingData;
    backend = nothing,
    concrete_core::Bool = false,
    group_repeated_algebraic_patterns = false
    )
    core = ExaModels.ExaCore(
        inf_model,
        data; 
        backend = backend, 
        concrete_core = concrete_core,
        group_repeated_algebraic_patterns = group_repeated_algebraic_patterns
    )
    return ExaModels.ExaModel(core)
end
function ExaModels.ExaModel(
    inf_model::InfiniteOpt.InfiniteModel;
    backend = nothing,
    group_repeated_algebraic_patterns = false,
    concrete_core::Bool = false
)
    return ExaModels.ExaModel(
        inf_model,
        ExaMappingData();
        backend = backend,
        group_repeated_algebraic_patterns = group_repeated_algebraic_patterns,
        concrete_core = concrete_core
    )
end
