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
    exaBackend = model.backend
    exaData = exaBackend.data
    exaModel = exaBackend.model

    # Test finite parameter mappings
    @test length(keys(exaBackend.data.finparam_mappings)) == 3
    @test InfiniteExaModels._check_mapping(x, exaBackend) === nothing
    @test InfiniteExaModels._check_mapping(y[1], exaBackend) === nothing
    @test InfiniteExaModels._check_mapping(y[2], exaBackend) === nothing
    xMapping = exaBackend.data.finparam_mappings[x]
    y1Mapping = exaBackend.data.finparam_mappings[y[1]]
    y2Mapping = exaBackend.data.finparam_mappings[y[2]]
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
    @variable(model, 0 ≤ v ≤ 100, Infinite(t))
    @parameter_function(model, pf == sin(t))
    @constraint(model, c1, v + pf ≤ 100)

    # Build the ExaModel backend
    build_transformation_backend!(model, model.backend)
    exaBackend = model.backend
    exaData = exaBackend.data
    exaModel = exaBackend.model

    # Test parameter function mappings
    @test !isempty(exaData.pfunc_info)
    @test exaData.pfunc_info[pf] isa Tuple{Symbol, ExaModels.Parameter}
    (alias, pfuncExa) = exaData.pfunc_info[pf]
    @test alias == :pf1
    @test pfuncExa isa ExaModels.Parameter
    @test pfuncExa.length == 5
    @test length(exaModel.θ) == 5
    @test exaModel.θ[pfuncExa.offset + 1 : pfuncExa.offset + pfuncExa.length] == sin.([0., 0.25, 0.5, 0.75, 1.0])

    # Test parameter function usage in constraints
    exaConstr = exaData.constraint_mappings[c1]
    @test exaConstr isa ExaModels.Constraint
    println(exaConstr)
end