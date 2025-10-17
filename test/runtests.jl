import Pkg
Pkg.add(url="https://github.com/infiniteopt/InfiniteOpt.jl/", rev="master")
using Test, InfiniteOpt, InfiniteExaModels, MadNLP, NLPModelsIpopt, ExaModels, Ipopt

println("************************************************")
println("                BEGINNING TESTS                 ")
println("************************************************\n")

@testset "Transcription Backend" begin include("transcription.jl") end
@testset "Solve Tests" begin include("solve.jl") end
@testset "InfiniteExaModelsMadNLP" begin include("madnlp.jl") end
@testset "InfiniteExaModelsIpopt" begin include("ipopt.jl") end

# TODO add more tests

println("\n************************************************")
println("                TESTING COMPLETE!                ")
println("************************************************")