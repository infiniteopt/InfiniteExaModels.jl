@testset "Finite Parameters" begin
    # Set up an InfiniteModel with finite parameters
    yVals = [20, 30]
    model = InfiniteModel(ExaTranscriptionBackend(NLPModelsIpopt.IpoptSolver))
    @infinite_parameter(model, t ∈ [0, 1], num_supports = 5)
    @finite_parameter(model, x == 42)
    @finite_parameter(model, y[i in 1:2] == yVals[i])
    @variable(model, 0 ≤ v ≤ 100, Infinite(t))

    # Build the ExaModel backend
    build_transformation_backend!(model, model.backend)
    exaBackend = InfiniteOpt.transformation_backend(model)
    @test exaBackend.core isa ExaModels.ExaCore
    exaModel = InfiniteOpt.transformation_model(model)
    exaData = InfiniteOpt.transformation_data(exaBackend)
    
    # Test finite parameter mappings
    @test length(keys(exaData.param_mappings)) == 3
    @test InfiniteExaModels._check_mapping(x, exaBackend) === nothing
    @test InfiniteExaModels._check_mapping(y[1], exaBackend) === nothing
    @test InfiniteExaModels._check_mapping(y[2], exaBackend) === nothing
    xMapping = exaData.param_mappings[x]
    y1Mapping = exaData.param_mappings[y[1]]
    y2Mapping = exaData.param_mappings[y[2]]
    @test xMapping isa ExaModels.Parameter{Tuple{Int}, Int}
    @test y1Mapping isa ExaModels.Parameter{Tuple{Int}, Int}
    @test y2Mapping isa ExaModels.Parameter{Tuple{Int}, Int}

    # Test finite parameter queries
    @test transformation_variable(x, exaBackend) == xMapping
    @test transformation_variable(y[1], exaBackend) == y1Mapping
    @test transformation_variable(y[2], exaBackend) == y2Mapping
    @test InfiniteExaModels._map_variable(x, x.index_type, 0, exaData) isa ExaModels.ParameterNode{Int}
    @test InfiniteExaModels._map_variable(y[1], y[1].index_type, 0, exaData) isa ExaModels.ParameterNode{Int}
    @test InfiniteExaModels._map_variable(y[2], y[2].index_type, 0, exaData) isa ExaModels.ParameterNode{Int}
    @test length(exaModel.θ) == 3
    @test exaModel.θ[xMapping.offset + 1] == 42
    @test exaModel.θ[y1Mapping.offset + 1] == yVals[1]
    @test exaModel.θ[y2Mapping.offset + 1] == yVals[2]
end

@testset "Parameter Functions" begin
    model = InfiniteModel(ExaTranscriptionBackend(NLPModelsIpopt.IpoptSolver))
    @infinite_parameter(model, t ∈ [0, 1], num_supports = 5)
    @infinite_parameter(model, s ∈ [2, 3], num_supports = 5)
    @variable(model, 0 ≤ v ≤ 100, Infinite(t))
    @parameter_function(model, pf == sin(t))
    @parameter_function(model, pf2 == (t, s) -> cos(t)*s)
    @constraint(model, c1, v + pf ≤ 100)
    @constraint(model, c2, v*2 + pf*pf2 ≤ 100)

    # Build the ExaModel backend
    build_transformation_backend!(model, model.backend)
    exaBackend = InfiniteOpt.transformation_backend(model)
    exaData = InfiniteOpt.transformation_data(exaBackend)
    exaModel = InfiniteOpt.transformation_model(model)

    # Test parameter function mappings (single infinite parameter)
    @test !isempty(exaData.param_mappings)
    @test exaData.param_mappings[pf] isa ExaModels.Parameter
    pfuncExa = exaData.param_mappings[pf]
    @test pfuncExa isa ExaModels.Parameter
    @test pfuncExa.length == 5
    expected = sin.([0., 0.25, 0.5, 0.75, 1.0])
    @test exaModel.θ[
        pfuncExa.offset + 1 : pfuncExa.offset + pfuncExa.length
        ] == expected 
    
    # test parameter function mappings (multiple infinite parameters)
    @test exaData.param_mappings[pf2] isa ExaModels.Parameter
    pfuncExa2 = exaData.param_mappings[pf2]
    @test pfuncExa2 isa ExaModels.Parameter
    @test pfuncExa2.length == 25
    @test length(exaModel.θ) == 30
    t_vals = [0., 0.25, 0.5, 0.75, 1.0]
    s_vals = [2, 2.25, 2.5, 2.75, 3]
    expected = vec([cos(t)*s for t in t_vals, s in s_vals])
    @test exaModel.θ[
        pfuncExa2.offset + 1 : pfuncExa2.offset + pfuncExa2.length
    ] == expected

    # Test parameter function usage in constraints
    exaConstr = exaData.constraint_mappings[c1]
    @test exaConstr isa ExaModels.Constraint

    exaConstr = exaData.constraint_mappings[c2]
    @test exaConstr isa ExaModels.Constraint
end

@testset "Objectives" begin
    model = InfiniteModel()
    @infinite_parameter(model, t in [0, 1], num_supports = 3)
    @infinite_parameter(model, x[1:2] in [-1, 1], independent = true, num_supports = 5)
    @variable(model, y, Infinite(t, x...))
    @variable(model, q, Infinite(t))
    @variable(model, z)
    x1_int = ∫(y^2, x[1])
    # Test ok objectives
    objs = [
        ∫(∫(x1_int, x[2]), t),
        ∫(∫(x1_int, x[2]) + 2q^2, t),
        ∫(∫(x1_int, x[2]) * sin(q), t),
        ∫(∫(x1_int, x[2]) + 2q, t),
        ∫(∫(x1_int, x[2]) + sin(q), t),
    ]
    for obj in objs
        @objective(model, Min, obj)
        @test ExaModel(model) isa ExaModel
    end
    # Test not so good objectives
    objs = [
        ∫(∫(x1_int^2, x[2]), t),
        ∫(∫(sin(x1_int), x[2]), t),
        ∫(∫(x1_int * x1_int^2.3 , x[2]), t),
        ∫(∫(sin(x1_int^2), x[2]), t),
        ∫(∫(x1_int + ∫(sin(y), x[1]), x[2]), t),
    ]
    for obj in objs
        @objective(model, Min, obj)
        @test_logs (:warn, InfiniteExaModels._ObjMeasureExpansionWarn) ExaModel(model) isa ExaModel
    end
end