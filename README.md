# Interfacing InfiniteOpt and ExaModels
The conversion API is in `InfiniteExaModels.jl` and the principal method is `exa_model` which creates
an `ExaModel` from an `InfiniteModel`. 

## Dependencies
The necessary packages are provided in `Project.toml`, but note that the compat hasn't been set up yet. 
Currently, the `master` branch of JuMP and the `jump_nlp` branch of InfiniteOpt are required.