# Excitonic BSE Model

The core of this package is to provide a general solution for establishing a excitonic Bethe-Salpeter model based on electronic Wannier functions, allowing users to customize the model.

## Construct

We provide four BSE models at present, they can all be constructed by a general function `BSE`:
```@docs
BSE
```

## The type of BSE Model

### :SU2

When any electronic state is doubly degenerated due to spin, excitonic states can be classified into the spin-singlet states and the spin-triplet states. In this case, you need to provide

|             |                                                                          |
| :---------- | :----------------------------------------------------------------------- |
| `sym`       | `:SU2`                                                                   |
| `TB`        | only contains one spin part, i.e. without spin                           |
| `Kernel`    | calculates $K^d$ and $K^x$ using electronic states that are without spin |
| `v` and `c` | the index of bands that are without spin                                 |

### :Sz

If we can separately provide the band structures of electrons with different spins, 
we can calculate the whole excitonic band structure based on the Sz of exciton. 
In this case, you need to provide

|                   |                                                    |
| :---------------- | :------------------------------------------------- |
| `sym`             | `:Sz`                                              |
| `TB_up`           | contains the spin-up electron structure            |
| `TB_up`           | contains the spin-down electron structure          |
| `Kernel`          | calculates $K^d$ and $K^x$ using electronic states |
| `v_up` and `c_up` | the index of bands of spin-up band                 |
| `v_dn` and `c_dn` | the index of bands of spin-down band               |

### :general

If any of the electronic states do not satisfy the spin degeneracy(noncollinear calculation), 
we have to calculate the whole excitonic band structure at the same time. 
In this case, you need to provide

|             |                                                                      |
| :---------- | :------------------------------------------------------------------- |
| `sym`       | `:general`                                                           |
| `TB`        | contains the whole electron structure, i.e. with spin                |
| `Kernel`    | calculates $K^d$ and $K^x$ using electronic states that contain spin |
| `v` and `c` | the index of bands that contain spin                                 |


## Kernel

Kernel is the most complex part of our BSE model, it is also the computational bottleneck because each element of the BSE Hamiltonian corresponds to a multiple summation. We provide Kernel implementations at [UJ approximation](kernel.md#U-approximation).

## Excitonic basis

In our package, every state is actually a numerical vector.
If you need to calculate some properties of a excitonic state, you may need know the basis order.
We provide two types of excitonic basis, $\left| vck \right\rangle$(for excitonic bloch function) and $\left| ijR \right\rangle$(for the periodic part of excitonic bloch function).

### $\left| vck \right\rangle$

If you already have a instance `bse` of type `AbstractBSE`, then you can get a `vckMap` by:

```julia
julia> bse.vckmap
vckMap(v = [1, 2], c = [3, 4], nk = 9)

julia> length(bse.vckmap) == 2 * 2 * 9
true

julia> size(bse.vckmap)
(2, 2, 9)

```

An object of `vckMap` supports bracket-based access to the basis's index.
```julia
julia> bse.vckmap[1]
(1, 3, 1)
```
This will give the `(vindex, cindex, kindex)`.
```julia
julia> bse.vckmap[1, 3, 1]
1
```
This is a reverse access compared to the previous one.

You may find you also need a actual kpoint, not only a `kindex`. 
The list of kpoints is `bse.kgrid`, you can get a kpoint by:

```julia
julia> bse.kgrid[kindex]
```

### $\left| ijR \right\rangle$

If you already have a instance `bse` of type `AbstractBSE`, then you can get a `ijRMap` by:

```julia
julia> bse.ijRmap
ijRMap(norb = 4, nR = 9)

julia> length(bse.ijRmap) == 4^2 * 9
true

julia> size(bse.ijRmap)
(4, 4, 9)

```

An object of `ijRMap` supports bracket-based access to the basis's index.
```julia
julia> bse.ijRmap[1]
(1, 1, 1)
```
This will give the `(i, j, Rindex)`.
```julia
julia> bse.ijRmap[1, 1, 1]
1
```
This is a reverse access compared to the previous one.

Also, you need a actual positive lattice vector, not only a `Rindex`. 
The list of positive lattice vectors is `bse.unitcell`, you can get a `R` by:
```julia
julia> bse.unitcell[Rindex]
```

## Hamiltonian

If you already have a instance `bse` of type `AbstractBSE`, then you can run:
```julia
julia> q = ReducedCoordinates(1, 2, 3)
julia> bse(q)
```
to get Hamiltonian of the BSE model.
For cluster, you can run:
```julia
julia> q = ReducedCoordinates(0, 0, 0)
julia> bse(q)
```
If the type of BSE is `:SU2`, `bse(q)` will return two matrixes at the same time:
```julia
julia> (Htriplet, Hsinglet) = bse(q);
```
If the type of BSE is `:Sz`, `bse(q)` will return two matrixes at the same time:
```julia
julia> (H0, Hp1, Hn1) = bse(q);
```
If the type of BSE is `:general`, `bse(q)` will return one matrix:
```julia
julia> H = bse(q);
```

It also support changing a given matrix directly:
```julia
julia> bse(Htriplet, Hsinglet, q) # For :SU2
```

## BAND

```@docs
BAND(kpoints::AbstractVector{<:ReducedCoordinates}, bse::WannierExcitonModel.BSESU2_wiTDA; kwargs...)
```

## Expectation of operator

You can also calculate the expectation value of an operator for each eigenstate, with the matrix representation of this operator.
This is the same as [tight-binding model](./tbmodel.md#Expectation-of-operator)

### spin operator

