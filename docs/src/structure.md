# Structure

We provide a group of types to help you build a crystalline structure.
We draw inspiration from [CrystallographyCore.jl](https://github.com/MineralsCloud/CrystallographyCore.jl) to design some types.
For most cases, these are only intermediate objects, and we use them to construct a *Model* type.

## Lattice

It can be used just like a matrix.

```@docs
Lattice(data::AbstractMatrix)
Lattice(𝐚::AbstractVector, 𝐛::AbstractVector, 𝐜::AbstractVector)
```

## Coordinates

They can be created and used just like a vector.

```@docs
ReducedCoordinates
CrystalCoordinates
CartesianCoordinates
```

**Example**

```@repl
using WannierExcitonModel

lattice = Lattice([
		   1.2 4.5 7.8
		   2.3 5.6 8.9
		   3.4 6.7 9.1
	   ])

a = ReducedCoordinates(1, 2, 3)
b = lattice * a
lattice \ b
```

## Cell

```@docs
Cell
```

### construct

```@docs
Cell(lattice, location::AbstractVector; kwargs...)
```

## ORBITAL

```@docs
ORBITAL
```

### construct

```@docs
ORBITAL(location::AbstractVector; kwargs...)
```
