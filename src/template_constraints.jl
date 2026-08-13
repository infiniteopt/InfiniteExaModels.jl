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
    return hash(-2, h), push!(refs, v), consts # -2 indicates a symbolic variable
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

function _exafy_iterated_expr(
    expr::Real,
    data::ExaMappingData
    )
    # TODO: finish
end
function _exafy_iterated_expr(
    expr::InfiniteOpt.GeneralVariableRef,
    data::ExaMappingData
    )
    # TODO: finish
end
function _exafy_iterated_expr(
    expr::Union{JuMP.GenericAffExpr{C, V}, JuMP.GenericQuadExpr{C, V}},
    data::ExaMappingData
    ) where {C, V}
    return _exafy_iterated_expr(convert(JuMP.GenericNonlinearExpr{V}, expr), data)
end
function _exafy_iterated_expr(
    expr::JuMP.GenericNonlinearExpr,
    data::ExaMappingData,
    
    )
    return _nl_op(expr.head)(_exafy_iterated_expr(a, data) for a in expr.args...)
end

# Given the 1st and last variable in an array of variables, create a new ExaModels.Variable 
# that represents the entire array
function _make_exa_variable(v1::ExaModels.Var, vf::ExaModels.Var, num_idxs::Int)
    num_vars = vf.i - v1.i + 1
    @assert num_vars == num_idxs
    return ExaModels.Variable((1:num_vars,), num_vars, v1.i, :x, nothing)
end
function _make_exa_variable(v1::ExaModels.Variable, vf::ExaModels.Variable, num_idxs::Int)
    @assert v1.size == vf.size
    num_vars = v1.length * num_idxs
    @assert vf.offset + vf.length - v1.offset == num_vars
    return ExaModels.Variable((v1.size..., 1:num_idxs), num_vars, v1.offset, :x, nothing)
end

function _process_array_object(
    obj::Vector{InfiniteOpt.GeneralVariableRef}, # TODO generalize to other arrays
    data::ExaMappingData
    )
    v_first = first(obj)
    v_last = last(obj)
    # TODO: what about infinite parameter arrays that are indexed across constraints?
    v_first.index_type in (InfiniteOpt.DependentParametersIndex, InfiniteOpt.IndependentParameterIndex) && return
    indexed_var = _make_exa_variable(data[v_first], data[v_last], length(obj))
    for v in obj
        @assert !haskey(data.indexed_var_to_exa_var, v)
        data.indexed_var_to_exa_var[v] = indexed_var
    end
    return
end
function _process_array_object(
    obj::Vector{InfiniteOpt.InfOptConstraintRef}, # TODO generalize to other arrays
    data::ExaMappingData
    )
    # encode each constraint
    h_ids = Vector{UInt}(undef, length(obj))
    vref_lists = Vector{Vector{InfiniteOpt.GeneralVariableRef}}(undef, length(obj))
    const_lists = Vector{Vector{Float64}}(undef, length(obj))
    for (i, cref) in enumerate(obj)
        con = JuMP.constraint_object(cref)
        expr = JuMP.jump_function(con)
        h_id, vrefs, consts = _encode_expr(expr)
        h_ids[i] = h_id
        vref_lists[i] = vrefs
        const_lists[i] = consts
    end
    # ensure we only have one unique constraint pattern in the array
    length(unique!(h_ids)) == 1 || return
    # determine which vrefs and consts change across the array
    vrefs1 = vref_lists[1]
    indexed_var_idxs = [any(l -> l[i] != vrefs1[i], vref_lists) for i in eachindex(vrefs1)]
    consts1 = const_lists[1]
    indexed_data_idxs = [any(l -> l[i] != consts1[i], const_lists) for i in eachindex(consts1)]
    # exafy the vrefs
    exafied_vrefs = Vector{Any}(undef, length(vrefs1))
    iterated_var_idx = 1
    for (i, vref) in enumerate(exafied_vrefs)
        if indexed_var_idxs[i]
            haskey(data.indexed_var_to_exa_var, vref) || return
            base_idxs = _index_params(vref, data)
            itr_alias = Symbol("itrvar$iterated_var_idx")
            var_idxs = (base_idxs..., ExaModels.DataSource()[itr_alias])
            exafied_vrefs[i] = data.indexed_var_to_exa_var[vref][var_idxs...]
            iterated_var_idx += 1
        else
            exafied_vrefs[i] = _exafy(vref, data)
        end
    end
    # exafy the consts
    exafied_consts = Vector{Any}(undef, length(consts1))

    # build the ExaModels graph for the constraint pattern

    # TODO: finish
    return
end
function _process_array_object(obj, data::ExaMappingData)
    return # fallback to ignore the unsupported/irrelevant object type
end

# Extract repeated constraint patterns using the object dictionary of the InfiniteModel
function _process_array_objects(
    core::ExaModels.ExaCore, 
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel
    )
    for (name, obj) in JuMP.object_dictionary(inf_model)
        _process_array_object(obj, data)
    end
end