# Interfacing InfiniteOpt and ExaModels
This provides an efficient bridge between InfiniteOpt and ExaModels. The principal method is `ExaTranscriptionBackend`
which creates an `AbstractTransformationBackend` for a `InfiniteModel`. 

## Dependencies
The necessary packages are provided in `Project.toml`. 
Currently, the `master` branch of InfiniteOpt is required.