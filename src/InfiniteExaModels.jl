module InfiniteExaModels

import InfiniteOpt, JuMP, ExaModels, SolverCore, NLPModels
import MathOptInterface as _MOI
import InfiniteOpt.TranscriptionOpt as _TO

include("infiniteopt_backend.jl")
include("operators.jl")
include("grouped_patterns.jl")
include("transform.jl")

export ExaMappingData, ExaTranscriptionBackend

end # end of module