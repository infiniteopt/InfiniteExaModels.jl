import Pkg
Pkg.add(url="https://github.com/infiniteopt/InfiniteOpt.jl/", rev="master")
using Test, InfiniteOpt, InfiniteExaModels, MadNLP, NLPModelsIpopt, ExaModels, Ipopt, Suppressor

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