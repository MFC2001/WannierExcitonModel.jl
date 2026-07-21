
"""
	KernelInterAction

This abstract type is used to calculate kernel term in excitonic BSE model.
You can design a concrete type of `KernelInterAction` according to your own needs.

For `YourKernel`, it should have specific fields:
- `nk::Int`: the number of kpoints;
- `kgrid::RedKgrid`: the kgrid used to calculate excitonic BSE model;
- `kgrid_Γ::RedKgrid`: with the same size of `kgrid` but without shift;
- `kgrid_addmap::Matrix{Int}`: satisfy `kgrid[i] + kgrid[j] = kgrid_Γ[kgrid_addmap[i, j]]`
- `kgrid_minusmap::Matrix{Int}`: satisfy `kgrid[i] - kgrid[j] = kgrid_Γ[kgrid_minusmap[i, j]]`

`YourKernel` should support these methods:

	(K::YourKernel)(::Val{:initialize})

It's used to calculate some interaction matrix at first, for example W(k′-k); it will be run when create a BSE model.

	(K::YourKernel)(q::ReducedCoordinates)

It's used to calculate interaction matrix related to q, for example V(q) or J(k+k′+q); it will be run when calculate BSE Hamiltonian.

	(K::YourKernel)(k′, k, ψ₁, ψ₂, ψ₃, ψ₄)

It's used to calculate kernel matrix elements ``K_{v'c'k',vck}``,
- `k′` and `k`: it's the index of kpoints in kgrid list;
- `ψ`: is a numerical vector; from left to right, in sequence: ``ψ_{c′}^{k′+q}``, ``ψ_{v′}^{k′}``, ``ψ_{v}^{k}``, ``ψ_{c}^{k+q}``.

!!! note
	You can change inputs, but use local operation carefully.
	Do not forget the minus of Kᵈ.

For the case only refer to q in a specified qgrid, i.e. `isqgrid` is true, this also means `k+q` is in the kgrid for any q needed.
In this case, `YourKernel` needs to support two more methods:

	(K::YourKernel)(::Val{:initialize_qgrid})

It's is a supplement to `(K::YourKernel)(::Val{:initialize})`; only run after `(K::YourKernel)(::Val{:initialize})` when create BSE model.

	(K::YourKernel)(kq_kindex::Vector{Int}, kΓq_kΓindex::Vector{Int}, q::ReducedCoordinates)

When `isqgrid` is true, BSEqgrid will not use (K::Kernel)(q::ReducedCoordinates) except when `q` equal to `[0, 0, 0]`, instead of using this function.
- `kq_kindex::Vector{Int}`: `kgrid[i] + q = kgrid[kq_kindex[i]]`;
- `kΓq_kΓindex::Vector{Int}`: `kgrid_Γ[i] + q = kgrid_Γ[kq_kindex[i]]`;
"""

abstract type KernelInterAction end

function (K::KernelInterAction)(sym::Symbol, args...; kwargs...)
	return K(Val(sym), args...; kwargs...)
end
function (K::KernelInterAction)(q::AbstractVector{<:Real})
	K.W(q)
	K.V(q)
	return nothing
end
function (K::KernelInterAction)(::Val{:buffer})
	buffer_Kᵈ = K.W(Val(:buffer))
	buffer_Kˣ = K.V(Val(:buffer))
	return buffer_Kᵈ, buffer_Kˣ
end
(K::KernelInterAction)(::Val{:buffer_Kᵈ}) = K.W(Val(:buffer))
(K::KernelInterAction)(::Val{:buffer_Kˣ}) = K.V(Val(:buffer))

include("./interaction.jl")
include("./_spinaware_Kernel.jl")
include("./_spinblind_Kernel.jl")
