using Pkg
Pkg.add(url="https://github.com/infiniteopt/InfiniteOpt.jl/", rev="master")
using Test, InfiniteOpt, InfiniteExaModels, MadNLP, NLPModelsIpopt, ExaModels

@test ExaTranscriptionBackend() isa ExaTranscriptionBackend

@testset "Finite Parameters" begin include("finite_parameters.jl") end

# TODO add tests