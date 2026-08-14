using Test, InfiniteOpt, InfiniteExaModels, MadNLP, NLPModelsIpopt
using ExaModels, Ipopt, Suppressor, NLPModels

# Global test functions
my_f(x) = x^2
my_∇f(x) = 2x
my_∇²f(x) = 2

println("************************************************")
println("                BEGINNING TESTS                 ")
println("************************************************\n")

@time @testset "Transcription Backend" begin include("transcription.jl") end
@time @testset "Solve Tests" begin include("solve.jl") end
@time @testset "InfiniteExaModelsMadNLP" begin include("madnlp.jl") end
@time @testset "InfiniteExaModelsIpopt" begin include("ipopt.jl") end

# TODO add more tests

println("\n************************************************")
println("                TESTING COMPLETE!                ")
println("************************************************")