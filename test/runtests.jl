import Pkg
Pkg.add(url="https://github.com/infiniteopt/InfiniteOpt.jl/", rev="master")
using Test, InfiniteOpt, InfiniteExaModels, MadNLP, NLPModelsIpopt, ExaModels, Ipopt

@testset "Transcription Backend" begin include("transcription.jl") end
@testset "Solve Tests" begin include("solve.jl") end

# TODO add more tests