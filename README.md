# Interfacing InfiniteOpt and ExaModels
The conversion API is in `InfiniteExaModels.jl` and the principal method is `exa_model` which creates
an `ExaModel` from an `InfiniteModel`. 

## Dependencies
The necessary packages are provided in `Project.toml`, but note that the compat is a little loose. 
Currently, the `jump_nlp` branch of InfiniteOpt is required.