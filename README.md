# WannierExcitonModel.jl

[![Build Status](https://github.com/MFC2001/WannierExcitonModel.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/MFC2001/WannierExcitonModel.jl/actions/workflows/CI.yml?query=branch%3Amaster)

[docs-url]: https://mfc2001.github.io/WannierExcitonModel.jl/

A Julia package solving excitonic Bethe-Salpeter equation(BSE) in electronic Wannier basis.

This package offers a general-purpose framework for constructing a BSE Hamiltonian in the electronic Wannier basis, 
with flexibility for user-defined model customization. 
In addition, it includes some post-processing utilities for further analysis and interpretation.

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
But if you purpose is to calculate the excitonic property, 
you only need to to clearly understand that the calculate process involves creating a `BSE` object and then using it.
See more details in Manual and Theory.

We require users to provide the Coulomb interaction matrix elements in Wannier basis.
which can be calculated by  `WannierInterAction.jl`.
