import Pkg
Pkg.add(url="https://github.com/infiniteopt/InfiniteOpt.jl/", rev="master")
using Test, InfiniteOpt, InfiniteExaModels, MadNLP, NLPModelsIpopt, ExaModels

@test ExaTranscriptionBackend() isa ExaTranscriptionBackend

@testset "Transcription Backend" begin include("transcription.jl") end
@testset "Parameter Function Methods" begin include("parameter_functions.jl") end

# TODO add tests