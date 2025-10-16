
export AbstractKernalInterAction
"""
	AbstractKernalInterAction

This abstract type is used to calculate kernal term in excitonic BSE model.
You can design a concrete type of `AbstractKernalInterAction` according to your own needs.

For `YourKernal`, it should have specific fields:
- `nk::Int`: the number of kpoints;
- `kgrid::RedKgrid`: the kgrid used to calculate excitonic BSE model;
- `kgrid_Γ::RedKgrid`: with the same size of `kgrid` but without shift;
- `kgrid_addmap::Matrix{Int}`: satisfy `kgrid[i] + kgrid[j] = kgrid_Γ[kgrid_addmap[i, j]]`
- `kgrid_minusmap::Matrix{Int}`: satisfy `kgrid[i] - kgrid[j] = kgrid_Γ[kgrid_minusmap[i, j]]`

`YourKernal` should support these methods:

	(K::YourKernal)(::Val{:initialize})

It's used to calculate some interaction matrix at first, for example W(k′-k); it will be run when create a BSE model.

	(K::YourKernal)(q::ReducedCoordinates)

It's used to calculate interaction matrix related to q, for example V(q) or J(k+k′+q); it will be run when calculate BSE Hamiltonian.

	(K::YourKernal)(k′, k, ψ₁, ψ₂, ψ₃, ψ₄)

It's used to calculate kernal matrix elements ``K_{v'c'k',vck}``,
- `k′` and `k`: it's the index of kpoints in kgrid list;
- `ψ`: is a numerical vector; from left to right, in sequence: ``ψ_{c′}^{k′+q}``, ``ψ_{v′}^{k′}``, ``ψ_{v}^{k}``, ``ψ_{c}^{k+q}``.

!!! note
	You can change inputs, but use local operation carefully.
	Do not forget the minus of Kᵈ.

For the case only refer to q in a specified qgrid, i.e. `isqgrid` is true, this also means `k+q` is in the kgrid for any q needed.
In this case, `YourKernal` needs to support two more methods:

	(K::YourKernal)(::Val{:initialize_qgrid})

It's is a supplement to `(K::YourKernal)(::Val{:initialize})`; only run after `(K::YourKernal)(::Val{:initialize})` when create BSE model.

	(K::YourKernal)(kq_kindex::Vector{Int}, kΓq_kΓindex::Vector{Int}, q::ReducedCoordinates)

When `isqgrid` is true, BSEqgrid will not use (K::Kernal)(q::ReducedCoordinates) except when `q` equal to `[0, 0, 0]`, instead of using this function.
- `kq_kindex::Vector{Int}`: `kgrid[i] + q = kgrid[kq_kindex[i]]`;
- `kΓq_kΓindex::Vector{Int}`: `kgrid_Γ[i] + q = kgrid_Γ[kq_kindex[i]]`;
"""
abstract type AbstractKernalInterAction end

function (K::AbstractKernalInterAction)(sym::Symbol, paras...)
	return K(Val(sym), paras...)
end


include("./UJ_trunction/UJ_trunction.jl")
include("./U_trunction/U_trunction.jl")
include("./cluster_UJ_trunction/cluster_UJ_trunction.jl")
