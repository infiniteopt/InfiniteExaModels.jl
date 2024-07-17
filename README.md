# InfiniteExaModels.jl
This package provides a transformation backend for [InfiniteOpt](https://github.com/infiniteopt/InfiniteOpt.jl)
such that InfiniteOpt models are efficiently transformed into [ExaModels](https://github.com/exanauts/ExaModels.jl)
via automated direct transcription. The underlying ExaModels models leverage recurrent algebraic structure to 
facilitate accelerated solution on CPUs and GPUs. Moreover, InfiniteOpt provides an intuitive interface that
automates transcription and drastically reduced model creation time relative to solving JuMP models via
ExaModels' `Optimizer` interface.

![Abstract](schematic.JPG)

## Installation
InfiniteExaModels is nearly ready for official release, but still needs some polishing touches first. In the meantime,
you can try it out by installing the developmental versions of InfiniteExaModels and InfiniteOpt 
(Julia `v1.9` or newer is required):
```julia
using Pkg
Pkg.add(url = "https://github.com/infiniteopt/InfiniteExaModels.jl", rev = "main")
Pkg.add(url = "https://github.com/infiniteopt/InfiniteOpt.jl", rev = "master")
```

## Usage
InfiniteExaModels primarily provides `ExaTranscriptionBackend` which can be passed to an `InfiniteModel` along
with a solver that is compliant with [JuliaSmoothOptimizers](https://github.com/JuliaSmoothOptimizers) standards.

### CPU Usage
Typical CPU workflows will use [Ipopt](https://github.com/JuliaSmoothOptimizers/NLPModelsIpopt.jl):
```julia
using InfiniteOpt, InfiniteExaModels, NLPModelsIpopt

model = InfiniteModel(ExaTranscriptionBackend(IpoptSolver))
```

### GPU Usage
Typical GPU workflows will use [MadNLP](https://github.com/MadNLP/MadNLP.jl), [CUDA](https://github.com/JuliaGPU/CUDA.jl]), 
and [CUDss](https://github.com/exanauts/CUDSS.jl) (a compatible Nvidia GPU is required):
```julia
using InfiniteOpt, InfiniteExaModels, MadNLP, CUDA # be sure to install CUDSS first as well

model = InfiniteModel(ExaTranscriptionBackend(MadNLPSolver, backend = CUDABackend()))
```

## Citation
If this is useful for your work please consider citing it:
```latex
@incollection{pulsipher2024scalable,
  title={Scalable Modeling of Infinite-Dimensional Nonlinear Programs with InfiniteExaModels.jl},
  author={Pulsipher, Joshua L and Shin, Sungho},
  booktitle={Computer Aided Chemical Engineering},
  volume={53},
  pages={3373--3378},
  year={2024},
  publisher={Elsevier}
}
```