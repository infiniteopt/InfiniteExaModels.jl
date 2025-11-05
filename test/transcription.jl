@testset "Mapping Initializers" begin
    # Setup model
    m = InfiniteModel()
    @infinite_parameter(m, t in [0, 1], num_supports = 5)
    @infinite_parameter(m, x in [-1, 1], num_supports = 3, 
                        derivative_method = OrthogonalCollocation(3))
    @variable(m, 1 >= y >= cos, Infinite(t))
    @variable(m, q == 42, Infinite(t, x))
    @variable(m, 2 <= w <= sin, Infinite(x), start = cos)
    y0 = y(0)
    y1 = y(1)
    set_start_value(y0, 0.5)
    delete_lower_bound(y1)
    set_upper_bound(y1, 0.8)
    q0 = q(0, x)
    q1 = q(1, x)
    set_start_value(q0, 10)
    fix(q1, 5)
    d1 = deriv(y, t)
    d2 = deriv(q, x, x)
    @variable(m, z, start = 10)
    c = ExaModels.ExaCore()
    data = ExaMappingData()
    # test building up supports
    @testset "_build_base_iterators" begin
        @test InfiniteExaModels._build_base_iterators(data, m) isa Nothing
        @test length(data.base_itrs) == 2
        @test sum(data.has_internal_supps) == 1
    end
    # test mapping finite variables
    @testset "_add_finite_variables" begin
        @test InfiniteExaModels._add_finite_variables(c, data, m) isa Nothing
        v = ExaModels.Var(1)
        @test data.finvar_mappings[z] == v
        @test c.x0[v.i] == 10
        @test c.lvar[v.i] == -Inf
        @test c.uvar[v.i] == Inf
    end
    # test mapping infinite variables
    @testset "_add_infinite_variables" begin
        @test InfiniteExaModels._add_infinite_variables(c, data, m) isa Nothing
        # test y variable mapping
        yvar = data.infvar_mappings[y]
        @test yvar.length == 5
        @test c.lvar[yvar.offset+1:yvar.offset+5] == cos.(range(0, 1, length=5))
        @test c.uvar[yvar.offset+1:yvar.offset+5] == ones(5)
        # test q variable mapping
        qvar = data.infvar_mappings[q]
        @test qvar.length == 25
        @test c.lvar[qvar.offset+1:qvar.offset+25] == fill(42, 25)
        @test c.uvar[qvar.offset+1:qvar.offset+25] == fill(42, 25)
        # test w variable mapping
        wvar = data.infvar_mappings[w]
        @test wvar.length == 5
        @test c.lvar[wvar.offset+1:wvar.offset+5] == fill(2, 5)
        @test c.uvar[wvar.offset+1:wvar.offset+5] == sin.(range(-1, 1, length=5))
        @test c.x0[wvar.offset+1:wvar.offset+5] == cos.(range(-1, 1, length=5))
        # test derivative mappings
        d1var = data.infvar_mappings[d1]
        @test d1var.length == 5
        @test data.infvar_mappings[d2].length == 25
        @test num_derivatives(m) == 3
    end
    # test semi-infinite variable mappings
    @testset "_add_semi_infinite_variables" begin
        @test InfiniteExaModels._add_semi_infinite_variables(c, data, m) isa Nothing
        @test length(data.semivar_info) == 2
        # test q0 variable mapping
        qvar = data.infvar_mappings[q]
        @test c.x0[qvar[1, 2].i] == 10
        # test q1 variable mapping
        @test c.lvar[qvar[5, 3].i] == 5
        @test c.uvar[qvar[5, 4].i] == 5
        @test c.x0[qvar[5, 2].i] == 0
    end
    # test point variable mappings
    @testset "_add_point_variables" begin
        @test InfiniteExaModels._add_point_variables(c, data, m) isa Nothing
        @test length(data.finvar_mappings) == 3
        # test y0 variable mapping
        yvar = data.finvar_mappings[y0]
        @test c.x0[yvar.i] == 0.5
        @test c.lvar[yvar.i] == cos(0)
        # test y1 variable mapping
        yvar = data.finvar_mappings[y1]
        @test c.lvar[yvar.i] == -Inf
        @test c.uvar[yvar.i] == 0.8
    end
end

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