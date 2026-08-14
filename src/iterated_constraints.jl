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

# Traverse expression in same order as _encode_expr and exafy it
function _exafy_iterated_expr(
    ::Real,
    vrefs::Vector{Any},
    consts::Vector{Any}
    )
    return popfirst!(consts)
end
function _exafy_iterated_expr(
    ::InfiniteOpt.GeneralVariableRef,
    vrefs::Vector{Any},
    consts::Vector{Any}
    )
    return popfirst!(vrefs)
end
function _exafy_iterated_expr(
    expr::Union{JuMP.GenericAffExpr{C, V}, JuMP.GenericQuadExpr{C, V}},
    vrefs::Vector{Any},
    consts::Vector{Any}
    ) where {C, V}
    return _exafy_iterated_expr(convert(JuMP.GenericNonlinearExpr{V}, expr), vrefs, consts)
end
function _exafy_iterated_expr(
    expr::JuMP.GenericNonlinearExpr,
    vrefs::Vector{Any},
    consts::Vector{Any}
    )
    return _nl_op(expr.head)((_exafy_iterated_expr(a, vrefs, consts) for a in expr.args)...)
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
function _make_exa_variable(v1::ExaModels.Parameter, vf::ExaModels.Parameter, num_idxs::Int)
    @assert v1.size == vf.size
    num_vars = v1.length * num_idxs
    @assert vf.offset + vf.length - v1.offset == num_vars
    return ExaModels.Parameter((v1.size..., 1:num_idxs), num_vars, v1.offset, nothing)
end

# get the repeated index of a variable based on its direct exaified variable reference
function _get_repeated_idx(em_var::ExaModels.Var, repeated_var::ExaModels.Variable)
    idx = em_var.i - repeated_var.offset + 1
    @assert idx in repeated_var.size[end]
    return idx
end
function _get_repeated_idx(em_var::ExaModels.Variable, repeated_var::ExaModels.Variable)
    idx = (em_var.offset - repeated_var.offset) ÷ em_var.length + 1
    @assert idx in repeated_var.size[end]
    return idx
end

function _process_variable_array(
    obj::Vector{InfiniteOpt.GeneralVariableRef}, # TODO generalize to other arrays
    data::ExaMappingData
    )
    v_first = first(obj)
    v_last = last(obj)
    # TODO: what about infinite parameter arrays that are indexed across constraints?
    v_first.index_type in (InfiniteOpt.DependentParametersIndex, InfiniteOpt.IndependentParameterIndex) && return false
    indexed_var = _make_exa_variable(data[v_first], data[v_last], length(obj))
    for v in obj
        @assert !haskey(data.indexed_var_to_exa_var, v)
        data.indexed_var_to_exa_var[v] = indexed_var
    end
    return true
end
function _process_variable_array(obj, data::ExaMappingData)
    return false
end

function _process_constraint_array(
    core::ExaModels.ExaCore,
    obj::Vector{InfiniteOpt.InfOptConstraintRef}, # TODO generalize to other arrays
    data::ExaMappingData
    )
    # encode each constraint
    h_ids = Vector{UInt}(undef, length(obj))
    vref_lists = Vector{Vector{InfiniteOpt.GeneralVariableRef}}(undef, length(obj))
    const_lists = Vector{Vector{Float64}}(undef, length(obj))
    sets = Vector{_MOI.AbstractSet}(undef, length(obj))
    for (i, cref) in enumerate(obj)
        isempty(JuMP.owner_model(cref).constraints[JuMP.index(cref)].measure_indices) || return false, core
        con = JuMP.constraint_object(cref)
        expr = JuMP.jump_function(con)
        h_id, vrefs, consts = _encode_expr(expr)
        h_ids[i] = h_id
        vref_lists[i] = vrefs
        const_lists[i] = consts
        sets[i] = JuMP.moi_set(con)
    end
    # ensure we only have one unique constraint pattern in the array
    length(unique!(h_ids)) == 1 || return false
    # determine which vrefs and consts change across the array
    vrefs1 = vref_lists[1]
    indexed_var_idxs = [any(l -> l[i] != vrefs1[i], vref_lists) for i in eachindex(vrefs1)]
    consts1 = const_lists[1]
    indexed_data_idxs = [any(l -> l[i] != consts1[i], const_lists) for i in eachindex(consts1)]
    # exafy the vrefs
    exafied_vrefs = Vector{Any}(undef, length(vrefs1))
    var_itr = Any[(;) for _ in 1:length(obj)]
    iterated_var_idx = 1
    for (i, vref) in enumerate(vrefs1)
        if indexed_var_idxs[i]
            haskey(data.indexed_var_to_exa_var, vref) || return false, core
            base_idxs = _index_params(vref, data)
            itr_alias = Symbol("itrvar$iterated_var_idx")
            var_idxs = (base_idxs..., ExaModels.DataSource()[itr_alias])
            src_var = data.indexed_var_to_exa_var[vref]
            exafied_vrefs[i] = src_var[var_idxs...]
            for j in 1:length(obj)
                em_var = data[vref_lists[j][i]]
                var_itr[j] = (; var_itr[j]..., itr_alias => _get_repeated_idx(em_var, src_var))
            end
            iterated_var_idx += 1
        else
            exafied_vrefs[i] = _exafy(vref, data)
        end
    end
    # exafy the consts
    exafied_consts = Vector{Any}(undef, length(consts1))
    const_itr = Any[(;) for _ in 1:length(obj)]
    iteratered_const_idx = 1
    for (i, c) in enumerate(consts1)
        if indexed_data_idxs[i]
            itr_alias = Symbol("itrconst$iteratered_const_idx")
            exafied_consts[i] = ExaModels.DataSource()[itr_alias]
            for j in 1:length(obj)
                const_itr[j] = (; const_itr[j]..., itr_alias => const_lists[j][i])
            end
            iteratered_const_idx += 1
        else
            exafied_consts[i] = c
        end
    end
    # build the ExaModels graph for the constraint pattern
    raw_expr = JuMP.jump_function(JuMP.constraint_object(first(obj)))
    em_expr = _finalize_expr(_exafy_iterated_expr(raw_expr, exafied_vrefs, exafied_consts))
    # process the iterator
    infinite_itr = _get_constraint_iterator(first(obj), data)
    finite_itr = [merge(var_itr[i], const_itr[i]) for i in 1:length(obj)]
    itr = vec([merge(i...) for i in Iterators.product(infinite_itr, finite_itr)])
    # add the constraints to the core
    lbs = Vector{Float64}(undef, length(obj))
    ubs = Vector{Float64}(undef, length(obj))
    for (i, s) in enumerate(sets)
        lbs[i], ubs[i] = _get_constr_bounds(s)
    end
    full_lbs = repeat(lbs, inner = length(infinite_itr))
    full_ubs = repeat(ubs, inner = length(infinite_itr))
    core, con = ExaModels.add_con(core, em_expr, itr, lcon = full_lbs, ucon = full_ubs)
    # save the constraint mappings
    inf_len = length(infinite_itr)
    for (i, cref) in enumerate(obj)
        base_idx = (i - 1) * inf_len + 1
        sliced_itr = itr[base_idx:base_idx + inf_len - 1]
        offset = con.offset + base_idx - 1
        data.constraint_mappings[cref] = ExaModels.Constraint(con.f, sliced_itr, offset, (1:inf_len,), nothing)
    end
    return true, core
end
function _process_constraint_array(core::ExaModels.ExaCore, obj, data::ExaMappingData)
    return false, core # fallback to ignore the unsupported/irrelevant object type
end

# Extract repeated constraint patterns using the object dictionary of the InfiniteModel
function _process_array_objects(
    core::ExaModels.ExaCore,
    data::ExaMappingData,
    inf_model::InfiniteOpt.InfiniteModel;
    print_info::Bool = true
    )
    # create ExaModels.Variable for each array of variables in the model
    # have to handle variables first before constraints since dictionary is unordered
    for (name, obj) in JuMP.object_dictionary(inf_model)
        success = _process_variable_array(obj, data)
        if success && print_info
            @info "Successfully processed finite indices for variable array: $name"
        end
    end
    # create the necessary metadata for each array of constraints in the model
    for (name, obj) in JuMP.object_dictionary(inf_model)
        success, core = _process_constraint_array(core, obj, data)
        if success && print_info
            @info "Successfully added finite indices for constraint array: $name"
        end
    end
    return core
end
