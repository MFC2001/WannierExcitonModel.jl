# WannierExcitonModel.jl

[![Build Status](https://github.com/MFC2001/WannierExcitonModel.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/MFC2001/WannierExcitonModel.jl/actions/workflows/CI.yml?query=branch%3Amaster)

[docs-url]: https://mfc2001.github.io/WannierExcitonModel.jl/

A Julia package solving Excitonic BSE problem based on Electronic Wannier model.

This package provides a general solution for establishing a BSE model based on electronic Wannier functions, allowing users to customize the model while also implementing several typical scenarios.

## Installation

The package can be installed with the Julia package manager.
From [the Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/), type `]` to enter
the [Pkg mode](https://docs.julialang.org/en/v1/stdlib/REPL/#Pkg-mode) and run:

```julia
pkg> add WannierExcitonModel
```

Or, equivalently, via [`Pkg.jl`](https://pkgdocs.julialang.org/v1/):

```julia
julia> import Pkg; Pkg.add("WannierExcitonModel")
```

## Documentation

-   [**DEVEL**][docs-url] - *documentation of the in-development version.*

## Quick start

This package provides many types and methods which may confuse you.
But if you purpose is to calculate the excitonic property, you only need to to clearly understand that the calculate process involves creating a `BSE` object and then using it.
See more details in Manual and Theory.

## In development

At present, we require users to provide the interaction between wannier bases.
For me, I only know [RESPACK](https://www.sciencedirect.com/science/article/pii/S001046552030391X) give this function, but it's a little bit unstable.
We are trying to calculate the interaction between wannier bases based on well-used open source software, such as [BerkeleyGW](https://berkeleygw.org/), by ourself.
