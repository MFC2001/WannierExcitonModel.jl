# Excitonic BSE Model

The core of this package is providing a general solution for establishing a excitonic Bethe-Salpeter model based on electronic Wannier functions, allowing users to customize the model while also implementing several typical scenarios (see more details in Theory).

For periodic case:

| Fields     | information                                                                               |
| :--------- | :---------------------------------------------------------------------------------------- |
| `TB`       | tight binding model                                                                       |
| `scissor`  | band gap correction                                                                       |
| `kgrid`    | list of kpoints `k`                                                                       |
| `unitcell` | list of positive lattice vector `R`                                                       |
| `vckmap`   | see behind                                                                                |
| `ijRmap`   | see behind                                                                                |
| `bandk`    | electronic band structure at each `k`, calculated when a BSE model is created             |
| `bandkq`   | electronic band structure at each `k+q`, calculated when calculated excitonic Hamiltonian |
| `Kernal`   | calculate $K^d$ and $K^x$ using electronic state which is without spin                    |

For cluster case:

| Fields    | information                                                                   |
| :-------- | :---------------------------------------------------------------------------- |
| `TB`      | tight binding model                                                           |
| `scissor` | band gap correction                                                           |
| `vcmap`   | see behind                                                                    |
| `ijmap`   | see behind                                                                    |
| `band`    | electronic band structure at each `k`, calculated when a BSE model is created |
| `Kernal`  | calculate $K^d$ and $K^x$ using electronic state which is without spin        |

## Construct

We provides four BSE models type at present, they all can be construct by a general function `BSE`:
```@docs
BSE
```

## The type of BSE Model

### :spinless

When any electronic state is doubly degenerate due to spin, excitonic states can be classified into the spin-singlet state and the spin-triplet state. In this case, you need to provide

|             |                                                                        |
| :---------- | :--------------------------------------------------------------------- |
| `sym`       | `:spinless`                                                            |
| `TB`        | only contain one spin part, i.e without spin                           |
| `Kernal`    | calculate $K^d$ and $K^x$ using electronic state which is without spin |
| `v` and `c` | the index of band which is without spin                                |

### :spinful

If any of the electronic states do not satisfy the spin degeneracy, we have to calculate the excitonic spin-singlet state and spin-triplet state at the same time. In this case, you need to provide

|             |                                                                      |
| :---------- | :------------------------------------------------------------------- |
| `sym`       | `:spinful`                                                           |
| `TB`        | contain the whole electron structure, i.e with spin                  |
| `Kernal`    | calculate $K^d$ and $K^x$ using electronic state which contains spin |
| `v` and `c` | the index of band which is with spin                                 |

### :cluster_spinless

This type is used when you model is a cluster and also spin degenerate. In this case, you only need to know setting `kgrid` is a invalid operation. Remem ber to set `period` of `TB` as `[0, 0, 0]`.

### :cluster_spinful

This type is used when you model is a cluster and also not spin degenerate. In this case, you only need to know setting `kgrid` is a invalid operation. Remem ber to set `period` of `TB` as `[0, 0, 0]`.

## Kernal

Kernal is the most complex part of our BSE model, it's also the computational bottlenack because each elements of the BSE Hamiltonian corresponds to a multiple summation. We provide two Kernel implementations under [UJ approximation](kernal.md#UJ-approximation) and the [U approximation](kernal.md#U-approximation) respectively. 
!!! note
    For spinful case, we require that the electronic wannier basis satisfy time-reversal symmetry, so that we can define which term is `U` and which term is `J`.

## Excitonic basis

In our package, every state is actually a numerical vector.
If you need to calculate some property of a excitonic state, you may need to know the basis order.
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

A object of `vckMap` supports bracket-based access to the basis's index.

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

You must find you also need a actual kpoint, not only a `kindex`. 
The list of kpoints is `bse.kgrid`, you can get a kpoint by:

```julia
julia> bse.kgrid[kindex]
```

For cluster, replace `vckmap` with `vcmap`:

```julia
julia> bse.vcmap
vcMap(v = [1, 2], c = [3, 4])

julia> length(bse.vcmap) == 2 * 2
true

julia> size(bse.vcmap)
(2, 2)

```

Also you can:

```julia
julia> bse.vcmap[1]
(1, 3)

julia> bse.vcmap[1, 3]
1

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

A object of `ijRMap` supports bracket-based access to the basis's index.

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
The list of positive lattice vectors is `bse.unitcell`, you can get a R by:

```julia
julia> bse.unitcell[Rindex]
```

For cluster, replace `ijRmap` with `ijmap`:
```julia
julia> bse.ijmap
ijMap(norb = 4)

julia> length(bse.ijmap) == 4^2
true

julia> size(bse.ijmap)
(4, 4)

```

Also you can:

```julia
julia> bse.ijmap[1]
(1, 1)
julia> bse.ijmap[1, 1]
1
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

If the type of BSE is `:spinless`, `bse(q)` will return two matrix at the same time:

```julia
julia> (Htriplet, Hsinglet) = bse(q);
```

If the type of BSE is `:spinful`, `bse(q)` will return one matrix:

```julia
julia> H = bse(q);
```

It also support changing a given matrix directly:

```julia
julia> bse(Htriplet, Hsinglet, q) # For :spinless
julia> bse(H, q) # For :spinful
```

## BAND

```@docs
BAND(kpoints::AbstractVector{<:ReducedCoordinates}, bse::WannierExcitonModel.BSEspinless; kwargs...)
```


## Expectation of operator

You can also calculate a operator's expectation value for every state, with the matrix representation of this operator.
This is the same as [tight-binding model](./tbmodel.md#Expectation-of-operator)

### spin operator

For the spinful BSE model, we provide some method to help you get the matrix representation of $\hat{S}_z$ and $\hat{S}^{2}$.

```julia
julia> numorb(TB)
6
julia> TB(:spinmat, [1, 2, 3]) # equal to TB(:spinmat, [1, 2, 3], [4, 5, 6])
```
its unit is $\frac{\hbar}{2}$.
