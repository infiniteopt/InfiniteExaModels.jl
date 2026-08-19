# Integer alias for types of InfiniteOpt modelling objects to use in hashing expressions
const _VariableTypeHashingInt = Dict(
    InfiniteOpt.FiniteParameterIndex => -2,
    InfiniteOpt.ParameterFunctionIndex => -3,
    InfiniteOpt.InfiniteVariableIndex => -4,
    InfiniteOpt.DerivativeIndex => -4,
    InfiniteOpt.SemiInfiniteVariableIndex => -4,
    InfiniteOpt.PointVariableIndex => -4,
    InfiniteOpt.FiniteVariableIndex => -5,
    InfiniteOpt.IndependentParameterIndex => -6,
    InfiniteOpt.DependentParameterIndex => -7,
    InfiniteOpt.MeasureIndex => -8,
)

## Extract the following from an expression:
# 1. A hash of the symbolic expression structure
# 2. A list of all variable references in the expression in the order they appear
# 3. A list of all constant values in the expression in the order they appear
function _encode_expr(expr::JuMP.AbstractJuMPScalar)
    return _encode_expr(expr, hash(:+), InfiniteOpt.GeneralVariableRef[], Float64[])
end
function _encode_expr(c::Real, h::UInt, refs, consts)
    return hash(-1, h), refs, push!(consts, c) # -1 indicates a symbolic constant
end
function _encode_expr(v::InfiniteOpt.GeneralVariableRef, h::UInt, refs, consts)
    return hash(_VariableTypeHashingInt[v.index_type], h), push!(refs, v), consts
end
function _encode_expr(
    expr::Union{JuMP.GenericAffExpr{C, V}, JuMP.GenericQuadExpr{C, V}},
    h::UInt,
    refs,
    consts
    ) where {C, V}
    return _encode_expr(convert(JuMP.GenericNonlinearExpr{V}, expr), h, refs, consts)
end
function _encode_expr(expr::JuMP.GenericNonlinearExpr, h::UInt, refs, consts) # TODO remove recursion
    h = hash((expr.head, length(expr.args)), h)
    for arg in expr.args
        h, _, _ = _encode_expr(arg, h, refs, consts)
    end
    return h, refs, consts
end

# Traverse expression in same order as _encode_expr and exafy it
function _exafy_grouped_expr(
    ::Real,
    vrefs::Vector{Any},
    consts::Vector{Any}
    )
    return popfirst!(consts)
end
function _exafy_grouped_expr(
    ::InfiniteOpt.GeneralVariableRef,
    vrefs::Vector{Any},
    consts::Vector{Any}
    )
    return popfirst!(vrefs)
end
function _exafy_grouped_expr(
    expr::Union{JuMP.GenericAffExpr{C, V}, JuMP.GenericQuadExpr{C, V}},
    vrefs::Vector{Any},
    consts::Vector{Any}
    ) where {C, V}
    return _exafy_grouped_expr(convert(JuMP.GenericNonlinearExpr{V}, expr), vrefs, consts)
end
function _exafy_grouped_expr(
    expr::JuMP.GenericNonlinearExpr,
    vrefs::Vector{Any},
    consts::Vector{Any}
    )
    return _nl_op(expr.head)((_exafy_grouped_expr(a, vrefs, consts) for a in expr.args)...)
end

# Print a message about a group
function _group_info_msg(group, msg)
    idxs = [JuMP.index(cref).value for cref in group]
    @info "$msg constraint group with indices: $(idxs)"
    return
end

# Get the grouped index of a variable based on its direct exaified variable reference
function _get_grouped_idx(em_var::ExaModels.Var, grouped_var::ExaModels.Variable)
    idx = em_var.i - grouped_var.offset
    @assert idx in grouped_var.size[end] && length(grouped_var.size) == 1
    return idx
end
function _get_grouped_idx(
    em_var::Union{ExaModels.Variable, ExaModels.Parameter},
    grouped_var::Union{ExaModels.Variable, ExaModels.Parameter}
    )
    idx = (em_var.offset - grouped_var.offset) ÷ em_var.length + 1
    @assert idx in grouped_var.size[end]
    return idx
end
function _get_grouped_idx(vref::InfiniteOpt.GeneralVariableRef, data::ExaMappingData)
    if vref.index_type in (InfiniteOpt.SemiInfiniteVariableIndex, InfiniteOpt.PointVariableIndex)
        vref = InfiniteOpt.infinite_variable_ref(vref)
    end
    em_var = data[vref]
    grouped_var = data.var_to_grouped_var[vref]
    return _get_grouped_idx(em_var, grouped_var)
end

# Given a candidate group of constraint, seek to merge together and add as a single constraint pattern to `core`
function _process_candidate_constraint_group(
    core::ExaModels.ExaCore,
    data::ExaMappingData,
    crefs::Vector{InfiniteOpt.InfOptConstraintRef},
    vref_lists::Vector{Vector{InfiniteOpt.GeneralVariableRef}},
    const_lists::Vector{Vector{Float64}},
    sets::Vector{_MOI.AbstractSet}
    )
    # determine which vrefs and consts change across the array
    vrefs1 = vref_lists[1]
    is_grouped_var = [any(l -> l[i] != vrefs1[i], vref_lists) for i in eachindex(vrefs1)]
    consts1 = const_lists[1]
    is_grouped_data = [any(l -> l[i] != consts1[i], const_lists) for i in eachindex(consts1)]
    # exafy the vrefs
    exafied_vrefs = Vector{Any}(undef, length(vrefs1))
    var_itr = Any[(;) for _ in 1:length(crefs)]
    group_var_idx = 1
    for (i, vref) in enumerate(vrefs1)
        if is_grouped_var[i]
            if !haskey(data.var_to_grouped_var, vref) 
                _group_info_msg(crefs, "Failure: $(vref) is not a grouped variable/parameter which prevents adding")
                return false, core
            end
            base_idxs = _index_params(vref, data)
            itr_alias = Symbol("grouped_var$group_var_idx")
            var_idxs = (base_idxs..., ExaModels.DataSource()[itr_alias])
            src_var = data.var_to_grouped_var[vref]
            exafied_vrefs[i] = src_var[var_idxs...]
            for j in 1:length(crefs)
                infvar = vref_lists[j][i]
                if data.var_to_grouped_var[infvar] != src_var
                    _group_info_msg(crefs, "Failure: Unable to correctly set up variable grouping for $(infvar) which prevents adding")
                    return false, core
                end
                var_itr[j] = (; var_itr[j]..., itr_alias => _get_grouped_idx(infvar, data))
            end
            group_var_idx += 1
        else
            exafied_vrefs[i] = _exafy(vref, data)
        end
    end
    # exafy the consts
    exafied_consts = Vector{Any}(undef, length(consts1))
    const_itr = Any[(;) for _ in 1:length(crefs)]
    grouped_const_idx = 1
    for (i, c) in enumerate(consts1)
        if is_grouped_data[i]
            itr_alias = Symbol("grouped_const$grouped_const_idx")
            exafied_consts[i] = ExaModels.DataSource()[itr_alias]
            for j in 1:length(crefs)
                const_itr[j] = (; const_itr[j]..., itr_alias => const_lists[j][i])
            end
            grouped_const_idx += 1
        else
            exafied_consts[i] = c
        end
    end
    # build the ExaModels graph for the constraint pattern
    raw_expr = JuMP.jump_function(JuMP.constraint_object(first(crefs)))
    em_expr = _finalize_expr(_exafy_grouped_expr(raw_expr, exafied_vrefs, exafied_consts))
    # process the iterator
    infinite_itr = _get_constraint_iterator(first(crefs), data)
    finite_itr = [merge(var_itr[i], const_itr[i]) for i in 1:length(crefs)]
    itr = vec([merge(i...) for i in Iterators.product(infinite_itr, finite_itr)])
    # add the constraints to the core
    lbs = Vector{Float64}(undef, length(crefs))
    ubs = Vector{Float64}(undef, length(crefs))
    for (i, s) in enumerate(sets)
        lbs[i], ubs[i] = _get_constr_bounds(s)
    end
    full_lbs = repeat(lbs, inner = length(infinite_itr))
    full_ubs = repeat(ubs, inner = length(infinite_itr))
    core, con = ExaModels.add_con(core, em_expr, itr, lcon = full_lbs, ucon = full_ubs)
    # save the constraint mappings
    inf_len = length(infinite_itr)
    for (i, cref) in enumerate(crefs)
        base_idx = (i - 1) * inf_len + 1
        sliced_itr = itr[base_idx:base_idx + inf_len - 1]
        offset = con.offset + base_idx - 1
        data.constraint_mappings[cref] = ExaModels.Constraint(con.f, sliced_itr, offset, (1:inf_len,), nothing)
    end
    return true, core
end

# Iterate over constraints in the InfiniteOpt model, group by algebraic pattern, and add to the ExaModels core
function _group_and_add_constraints(
    core::ExaModels.ExaCore,
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    # set up dictionaries for tracking patterns
    hash_to_patterns = Dict{UInt, Tuple{Vector{Vector{InfiniteOpt.GeneralVariableRef}}, Vector{Vector{Float64}}, Vector{_MOI.AbstractSet}}}()
    hash_to_constrs = Dict{UInt, Vector{InfiniteOpt.InfOptConstraintRef}}()
    # iterate over constraints and group by hashed algebraic pattern
    for cref in JuMP.all_constraints(inf_model)
        InfiniteOpt.is_variable_domain_constraint(cref) && continue
        isempty(JuMP.owner_model(cref).constraints[JuMP.index(cref)].measure_indices) || continue # TODO: temporary restriction
        expr = JuMP.jump_function(JuMP.constraint_object(cref))
        expr isa JuMP.AbstractJuMPScalar || continue
        h, vrefs, consts = _encode_expr(expr)
        if haskey(hash_to_patterns, h)
            push!(hash_to_patterns[h][1], vrefs)
            push!(hash_to_patterns[h][2], consts)
            push!(hash_to_patterns[h][3], JuMP.moi_set(JuMP.constraint_object(cref)))
            push!(hash_to_constrs[h], cref)
        else
            hash_to_patterns[h] = ([vrefs], [consts], [JuMP.moi_set(JuMP.constraint_object(cref))])
            hash_to_constrs[h] = [cref]
        end
    end
    # process each grouped pattern (requiring at least 2 constraints to be grouped)
    for (h, crefs) in hash_to_constrs
        if length(crefs) < 2
            continue
        end
        success = _process_candidate_constraint_group(core, data, crefs, hash_to_patterns[h]...)
        if success
            _group_info_msg(crefs, "Successfully added")
        end
    end
    return core
end