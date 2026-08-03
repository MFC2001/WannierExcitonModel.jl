# Tight Binding Model

This package provides a group of functions to process tight-binding model.
The core object is `TightBindModel`, 

```@docs
TightBindModel
```

## Construct

You can construct a `TightBindModel` by

```@docs
TightBindModel(;kwargs...)
```

We emphasize that there are three parts of `TightBindModel` are the important. They are:

|                         | common file format      |
| :---------------------- | :---------------------- |
| lattice basis vector    | `POSCAR.vasp`           |
| orbital center location | `wannier90_centres.xyz` |
| hopping terms           | `wannier90_hr.dat`      |

The simplest method to construct a `TightBindModel` is to run:

```julia
julia> TightBindModel(;	hops = "wannier90_hr.dat", cell = "POSCAR.vasp", orbital = "wannier90_centres.xyz", period)
```

Remember to set `period`.

## Hamiltonian

If you already have a instance `TB` of type `TightBindModel`, then you can run:
```julia
julia> k = ReducedCoordinates(1, 2, 3)
julia> TB(k)
```
to get Hamiltonian of the tight-binding model. Or you may want to store Hamiltonian matrix to a existing variable:
```julia
julia> k = ReducedCoordinates(1, 2, 3)
julia> H = Matrix{Float64}(undef, numorb(TB), numorb(TB))
julia> TB(H, k)
```

When you need the Hamiltonian under atom guage, whose eigen vectors represent $\left| u_{n\mathbf{k}} \right\rangle$(see more details in Theory), you can run:

```julia
julia> TB(k, TB.orb_location)
julia> TB(H, k, TB.orb_location)
```

Under atom guage, you can also get the partial derivative of Hamiltonian from 

```julia
julia> TB(:partial, k)
julia> TB(:partial, H, k)
julia> size(H) == (numorb(TB), numorb(TB), 3)
true
```

## BAND

The function `BAND` is the core of this package.

```@docs
BAND(kpoints::AbstractVector{<:ReducedCoordinates}, TB::TightBindModel; kwargs...)
```
```@docs
ExtractU(kpoints, band, TB::TightBindModel)
```

## Expectation of operator

You can calculate the expectation value of an operator for each eigenstate, with the matrix representation of this operator.

```@docs
ExpectationValue(operator::AbstractMatrix{<:Number}, band::AbstractVector{<:WannierExcitonModel.Eigen})
```

### spin operator

For the simplest case where the tight-binding model contains spin, we can obtain the $\hat{S}_{z}$ matrix by:

```julia
julia> numorb(TB)
6
julia> TB(:spinmat, [1, 2, 3]) # equal to TB(:spinmat, [1, 2, 3], [4, 5, 6])
```

The unit is $\frac{\hbar}{2}$.