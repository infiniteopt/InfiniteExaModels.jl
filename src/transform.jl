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
    info::JuMP.VariableInfo{<:Real, <:Real, <:Real, <:Real},
    itrs = nothing
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
    lb = fill(-Inf, dims)
    ub = fill(Inf, dims)
    start = fill(0.0, dims)
    # iterate over all support combinations and fill in the arrays
    for i in Iterators.product(itrs...)
        supp = [s for nt in i for s in Iterators.drop(values(nt), 1)]
        if info.has_fix
            val = _process_value(info.fixed_value, supp)
            lb[first.(i)...] = val
            ub[first.(i)...] = val
        end
        if info.has_lb
            lb[first.(i)...] = _process_value(info.lower_bound, supp)
        end
        if info.has_ub
            ub[first.(i)...] = _process_value(info.upper_bound, supp)
        end
        if info.has_start
            start[first.(i)...] = _process_value(info.start, supp)
        end
    end
    return lb, ub, start
end

# Add all the finite variables from an InfiniteModel to a ExaCore
function _add_finite_variables(
    core::ExaModels.ExaCore, 
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for vref in JuMP.all_variables(inf_model, InfiniteOpt.FiniteVariable)
        info = InfiniteOpt.core_object(vref).info # JuMP.VariableInfo
        _ensure_continuous(info)
        lb, ub, start = _get_variable_bounds_and_start(info)
        new_var = ExaModels.variable(core, 1, start = start, lvar = lb, uvar = ub)
        data.finvar_mappings[vref] = new_var[1]
    end
    return
end

# Add all the finite parameters from an InfiniteModel to a ExaCore
function _add_finite_parameters(
    core::ExaModels.ExaCore, 
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for pref in JuMP.all_variables(inf_model, InfiniteOpt.FiniteParameter)
        param_val = InfiniteOpt.parameter_value(pref)
        data.param_mappings[pref] = ExaModels.parameter(core, [param_val])
    end
    return
end

# Add all the infinite variables (and derivatives) from an InfiniteModel to a ExaCore
function _add_infinite_variables(
    core::ExaModels.ExaCore, 
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    # Get the raw variables
    ivrefs = JuMP.all_variables(inf_model, InfiniteOpt.InfiniteVariable)
    drefs = InfiniteOpt.all_derivatives(inf_model)
    # process lower order derivatives first if needed
    for dref in drefs
        var = InfiniteOpt.core_object(dref)
        method = InfiniteOpt.derivative_method(dref)
        if !InfiniteOpt.allows_high_order_derivatives(method) && var.order > 1
            for o in var.order-1:-1:1
                if !haskey(inf_model.deriv_lookup, (var.variable_ref, var.parameter_ref, o))
                    info = InfiniteOpt._DefaultDerivativeInfo
                    new_d = InfiniteOpt.Derivative(info, var.variable_ref, var.parameter_ref, o)
                    # TODO: avoid changing the InfiniteModel
                    new_dref = InfiniteOpt.add_derivative(inf_model, new_d)
                    push!(ivrefs, new_dref)
                end
            end
        end
    end
    # now process and add each infinite variable
    for vref in append!(ivrefs, drefs)
        # retrieve basic information
        info = InfiniteOpt.core_object(vref).info # JuMP.VariableInfo
        _ensure_continuous(info)
        group_idxs = InfiniteOpt.parameter_group_int_indices(vref)
        # prepare the bounds and start values
        itrs = map(i -> data.base_itrs[i], group_idxs)
        lb, ub, start = _get_variable_bounds_and_start(info, itrs)
        # create the ExaModels variable
        dims = Tuple(length(itr) for itr in itrs)
        new_var = ExaModels.variable(core, dims...; start = start, lvar = lb, uvar = ub)
        data.infvar_mappings[vref] = new_var
    end
    return 
end

# Process all the parameter function from an InfiniteModel and add to an ExaCore
function _add_parameter_functions(
    core::ExaModels.ExaCore,
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )  
    for pfref in InfiniteOpt.all_parameter_functions(inf_model)
        # gather the basic information
        group_idxs = InfiniteOpt.parameter_group_int_indices(pfref)
        pfunc = InfiniteOpt.core_object(pfref)
        # compute the value for each support combination and store
        itrs = map(i -> data.base_itrs[i], group_idxs)
        dims = Tuple(length(itr) for itr in itrs)
        vals = Array{Float64}(undef, dims...)
        for i in Iterators.product(itrs...)
            supp = [s for nt in i for s in Iterators.drop(values(nt), 1)]
            vals[first.(i)...] = pfunc(supp)
        end
        # Register the parameter function values in the ExaCore & mapping data
        data.param_mappings[pfref] = ExaModels.parameter(core, vals)
    end
    return
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
    return data.infvar_mappings[ivref][idxs...]
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
        data.finvar_mappings[vref] = pt
        # update the bounds and start value if needed
        info = InfiniteOpt.core_object(vref).info # InfiniteOpt.RestrictedDomainInfo
        _update_bounds_and_start(core, info, pt)
    end
    return 
end

# Map variable references based on their underlying type (used by `_exafy`)
function _map_variable(vref, ::Type{InfiniteOpt.FiniteVariableIndex}, par_src, data)
    return data.finvar_mappings[vref]
end
function _map_variable(vref, ::Type{InfiniteOpt.PointVariableIndex}, par_src, data)
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
    par_src,
    data
    ) where V <: Union{InfiniteOpt.InfiniteVariableIndex, InfiniteOpt.DerivativeIndex}
    group_idxs = InfiniteOpt.parameter_group_int_indices(vref)
    idx_pars = (par_src[data.group_alias[i]] for i in group_idxs)
    return data.infvar_mappings[vref][idx_pars...]
end
function _map_variable(vref, ::Type{InfiniteOpt.SemiInfiniteVariableIndex}, par_src, data)
    if !haskey(data.semivar_info, vref)
        _process_semi_infinite_var(vref, data)
    end
    ivar, inds = data.semivar_info[vref]
    idx_pars = (i isa Int ? i : par_src[i] for i in inds)
    return ivar[idx_pars...]
end
function _map_variable(vref, ::Type{<:InfiniteOpt.InfiniteParameterIndex}, par_src, data)
    return par_src[data.param_alias[vref]]
end
function _map_variable(vref, ::Type{InfiniteOpt.FiniteParameterIndex}, par_src, data)
    return data.param_mappings[vref][1]
end
function _map_variable(vref, ::Type{InfiniteOpt.ParameterFunctionIndex}, par_src, data)
    group_idxs = InfiniteOpt.parameter_group_int_indices(vref)
    idx_pars = (par_src[data.group_alias[i]] for i in group_idxs)
    return data.param_mappings[vref][idx_pars...]
end
function _map_variable(vref, IdxType, par_src, data)
    error("Unable to add `$vref` to an ExaModel, it's index type `$IdxType`" *
          " is not yet supported by InfiniteExaModels.")
end

# Convert as InfiniteOpt expression into a ExaModel expression using the ParIndexed `par_src`
function _exafy(vref::InfiniteOpt.GeneralVariableRef, par_src, data)
    return _map_variable(vref, vref.index_type, par_src, data)
end
function _exafy(c::Number, par_src, data)
    return c
end
function _exafy(
    aff::JuMP.GenericAffExpr{C, InfiniteOpt.GeneralVariableRef}, 
    par_src, 
    data
    ) where {C}
    c = JuMP.constant(aff)
    if !isempty(aff.terms)
        ex = sum(begin
            v_ex = _exafy(v, par_src, data)
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
    par_src, 
    data
    ) where {C}
    aff = _exafy(quad.aff, par_src, data)
    if !isempty(quad.terms)
        ex = sum(begin 
            if v1 == v2
                v_ex = _exafy(v1, par_src, data) 
                isone(c) ? abs2(v_ex) : c * abs2(v_ex)
            else
                v1_ex = _exafy(v1, par_src, data) 
                v2_ex = _exafy(v2, par_src, data) 
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
    par_src, 
    data
    )
    return _nl_op(nl.head)((_exafy(a, par_src, data) for a in nl.args)...)
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
    error("Constraint set `$set` is not supported by InfiniteExaModels.")
end

# Check if NamedTuple iterator respects the restriction
function _support_in_restriction(restriction, itr, data)
    supp = [itr[data.param_alias[p]] for p in restriction.parameter_refs]
    return restriction(supp)
end

# Add all the constraints from an InfiniteModel to an ExaCore
function _add_constraints(
    core::ExaModels.ExaCore, 
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for cref in JuMP.all_constraints(inf_model)
        # skip if the constraint is a variable bound or type
        InfiniteOpt.is_variable_domain_constraint(cref) && continue
        # parse the basic information
        constr = JuMP.constraint_object(cref)
        if isempty(inf_model.constraints[JuMP.index(cref)].measure_indices)
            expr = JuMP.jump_function(constr)
        else
            @warn "Constrained measures can lead to poor performance with ExaModels."
            expr = InfiniteOpt.expand_measures(JuMP.jump_function(constr), inf_model)
        end
        set = JuMP.moi_set(constr)
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
            filter!(i -> _support_in_restriction(restriction, i, data), itr)
        end
        # create the ExaModels expression tree based on expr
        par_src = ExaModels.ParSource()
        em_expr = _finalize_expr(_exafy(expr, par_src, data))
        # get the constraint bounds
        lb, ub = _get_constr_bounds(set)
        # create the ExaModels constraint
        con = ExaModels.constraint(core, em_expr, itr, lcon = lb, ucon = ub)
        data.constraint_mappings[cref] = con
    end
    return
end

# Make dispatch type to pass the data needed by `make_reduced_expr`
struct _DerivReductionBackendInfo <: InfiniteOpt.AbstractTransformationBackend
    data::ExaMappingData
    par_src::ExaModels.ParSource
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
    par_src = dispatch_data.par_src
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
                par_src[i] 
            end
            end for i in inds)
        return ivar[idx_pars...]
    else # either an infinite variable or a derivative variable
        group_idxs = InfiniteOpt.parameter_group_int_indices(vref)
        idx_pars = (begin 
            g_alias = data.group_alias[i]
            if g_alias == alias
                idx
            else
                par_src[g_alias]
            end 
            end for i in group_idxs)
        return data.infvar_mappings[vref][idx_pars...]
    end
    return
end

# Add the approximation equations for each derivative variable
function _add_derivative_approximations(
    core::ExaModels.ExaCore, 
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for dref in InfiniteOpt.all_derivatives(inf_model)
        # gather the derivative information
        vref = InfiniteOpt.derivative_argument(dref) 
        pref = InfiniteOpt.operator_parameter(dref)
        order = InfiniteOpt.derivative_order(dref)
        method = InfiniteOpt.derivative_method(dref)
        # take care of derivative nesting if appropriate
        if !InfiniteOpt.allows_high_order_derivatives(method) && order > 1
            d_idx = inf_model.deriv_lookup[vref, pref, order - 1]
            vref = InfiniteOpt.GeneralVariableRef(inf_model, d_idx)
        end
        # gather the needed infinite parameter data
        group_idxs = InfiniteOpt.parameter_group_int_indices(vref)
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
        idxs, arg_itrs... = InfiniteOpt.derivative_expr_data(dref, order, supps, method)
        # make the iterator
        aliases = Tuple(Symbol("d_arg$i") for i in eachindex(arg_itrs))
        pref_itr = [(; srt_itr[i]..., zip(aliases, args)...) for (i, args...) in zip(idxs, arg_itrs...)]
        if length(group_idxs) > 1
            itrs = [g == pref_group ? pref_itr : data.base_itrs[g] for g in group_idxs]
            itr = [merge(i...) for i in Iterators.product(itrs...)]
        else
            itr = pref_itr
        end
        # make the ExaModel expression tree
        par_src = ExaModels.ParSource()
        em_expr = InfiniteOpt.make_indexed_derivative_expr(
            dref, 
            vref, 
            pref, 
            order, 
            par_src[data.group_alias[pref_group]], 
            supps, 
            _DerivReductionBackendInfo(data, par_src),
            method,
            (par_src[a] for a in aliases)...
            )
        # add the constraint
        ExaModels.constraint(core, em_expr, itr)
    end
    return
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
        # make the constraints for each infinite variable
        for vidx in vidxs
            vref = InfiniteOpt.GeneralVariableRef(inf_model, vidx)
            group_idxs = InfiniteOpt.parameter_group_int_indices(vref)
            aliases = (data.group_alias[g] for g in group_idxs)
            itrs = (g == pref_group ? pref_itr : data.base_itrs[g] for g in group_idxs)
            itr = vec([merge(i...) for i in Iterators.product(itrs...)])
            par_src = ExaModels.ParSource()
            idx_pars1 = (a == pref_alias ? par_src[:i1] : par_src[a] for a in aliases)
            idx_pars2 = (a == pref_alias ? par_src[:i2] : par_src[a] for a in aliases)
            ivar = data.infvar_mappings[vref]
            em_expr = ivar[idx_pars1...] - ivar[idx_pars2...]
            ExaModels.constraint(core, em_expr, itr)
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
    em_expr = _finalize_expr(_exafy(expr, (;), data))
    ExaModels.objective(core, em_expr, [(;)])
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
function _add_objective_aff_term(core, coef, vref, data)
    return _add_objective_aff_term(core, coef, vref, vref.index_type, data)
end
function _add_objective_aff_term(core, coef, vref, ::Type{InfiniteOpt.MeasureIndex}, data)
    # process the measure structure recursively as needed
    mexpr, itr = _process_measure_sum(vref, data)
    # prepare the examodel expression tree
    par_src = ExaModels.ParSource()
    em_expr = par_src.c * _exafy(coef * mexpr, par_src, data)
    # add the term to the objective
    ExaModels.objective(core, _finalize_expr(em_expr), itr)
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
    data::ExaMappingData, 
    inf_model::InfiniteOpt.InfiniteModel
    )
    vrefs = InfiniteOpt.all_expression_variables(expr)
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
    data::ExaMappingData, 
    ::InfiniteOpt.InfiniteModel
    )
    _add_objective_aff_term(core, 1.0, vref, data)
    return
end
function _add_objective(
    core::ExaModels.ExaCore,
    aff::JuMP.GenericAffExpr,
    data::ExaMappingData,
    ::InfiniteOpt.InfiniteModel
    )
    # TODO should we check if there are a lot of terms?
    for (coef, vref) in JuMP.linear_terms(aff)
        _add_objective_aff_term(core, coef, vref, data)
    end
    c = JuMP.constant(aff)
    if !iszero(c)
        ExaModels.objective(core, ExaModels.Null(c))
    end
    return
end
function _add_objective(
    core::ExaModels.ExaCore,
    quad::InfiniteOpt.GenericQuadExpr,
    data::ExaMappingData, 
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

# Fill an ExaCore and ExaMappingData with an InfiniteModel
function build_exa_core!(
    core::ExaModels.ExaCore, 
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    # initial setup
    _build_base_iterators(data, inf_model)
    # add the variables and appropriate mappings
    _add_finite_parameters(core, data, inf_model)
    _add_finite_variables(core, data, inf_model)
    _add_infinite_variables(core, data, inf_model) # includes derivatives
    _add_parameter_functions(core, data, inf_model)
    _add_semi_infinite_variables(core, data, inf_model)
    _add_point_variables(core, data, inf_model)
    # add the constraints
    _add_constraints(core, data, inf_model)
    _add_derivative_approximations(core, data, inf_model)
    _add_collocation_restrictions(core, data, inf_model)
    # add the objective if there is one
    expr = JuMP.objective_function(inf_model)
    sense = JuMP.objective_sense(inf_model)
    if sense != _MOI.FEASIBILITY_SENSE
        _add_objective(core, expr, data, inf_model)
    end
    return
end

"""

"""
function ExaModels.ExaCore(
    inf_model::InfiniteOpt.InfiniteModel,
    data::ExaMappingData;
    backend = nothing
    )
    # TODO add support for other float types once InfiniteOpt does
    minimize = JuMP.objective_sense(inf_model) == _MOI.MIN_SENSE
    core = ExaModels.ExaCore(; backend = backend, minimize = minimize)
    build_exa_core!(core, data, inf_model)
    return core
end

"""

"""
function ExaModels.ExaModel(
    inf_model::InfiniteOpt.InfiniteModel,
    data::ExaMappingData;
    backend = nothing
    )
    core = ExaModels.ExaCore(inf_model, data; backend = backend)
    return ExaModels.ExaModel(core)
end
function ExaModels.ExaModel(inf_model::InfiniteOpt.InfiniteModel; backend = nothing)
    return ExaModels.ExaModel(inf_model, ExaMappingData(), backend = backend)
end
