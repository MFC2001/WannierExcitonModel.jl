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

The workflow has four steps; each one is demonstrated by a runnable script in
[`examples/`](examples), which is executed from its own directory containing the
`unitcell/` input files.
The following is a simple template.

```julia
using WannierExcitonModel, JLD2

# 1. tight-binding model from wannier90 output
TB = TightBindModel(;
    hops = "unitcell/wannier90_hr.dat",
    cell = "unitcell/POSCAR",
    orbital = "unitcell/wannier90_centres.xyz",
    period = [1, 1, 0],
)

# 2. screened and bare electron-electron interaction
Kᵈ = PairInterAction(load("unitcell/U_1_mirror.jld2")["Uwithlr"],
                     load("unitcell/pair_1.jld2")["interaction"])
Kˣ = PairInterAction(load("unitcell/U_bare_mirror.jld2")["Uwithlr"],
                     load("unitcell/pair_bare.jld2")["interaction"])

# 3. BSE Hamiltonian on a k-grid
bse = BSE(MonkhorstPack([36, 36, 1]), TB, Kᵈ, Kˣ; v = 13:14, c = 15:16, scissor = 1.2)

# 4. exciton band structure along a k-path
kline = Kline(30, bse)
band = BAND(kline, bse)
```

Spin-polarized systems are handled by the `BSE(kgrid, TB_up, TB_dn, ...; v_up, c_up, v_dn, c_dn)`
variant, and noncollinear ones by passing `channel = 1`. See
[`examples/MoS2/bseband.jl`](examples/MoS2/bseband.jl),
[`examples/LiF/excband.jl`](examples/LiF/excband.jl) and
[`examples/Graphone/bseband.jl`](examples/Graphone/bseband.jl) for complete cases;
post-processing such as Berry curvature and its integration is shown in
[`examples/eleBC.jl`](examples/eleBC.jl) and [`examples/MoS2/excBC.jl`](examples/MoS2/excBC.jl).
