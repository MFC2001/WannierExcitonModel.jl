
export AbstractKernalInterAction
"""
	AbstractKernalInterAction

	struct Kernal <: AbstractKernalInterAction
		nk::Int
		kgrid::RedKgrid
		kgrid_Γ::RedKgrid
		kgrid_addmap::Matrix{Int}
		kgrid_minusmap::Matrix{Int}
		...
	end
	kgrid[i] + kgrid[j] = kgrid_Γ[kgrid_addmap[i, j]]
	kgrid[i] - kgrid[j] = kgrid_Γ[kgrid_minusmap[i, j]]
	kgrid is your basis, kgrid_Γ don't have shift with the same size of kgrid.

	You need to overload these three functions for your `Kernal`.
	function (K::Kernal)(::Val{:initialize}) end
	it's used to calculate some interaction matrix at first, for example W(k′-k).
	function (K::Kernal)(q::ReducedCoordinates) end
	it's used to calculate interaction matrix related to q, for example V(q) or J(k+k′+q)
	function (K::Kernal)(k′, k, ψ_{c′}^{k′+q}, ψ_{v′}^{k′}, ψ_{v}^{k}, ψ_{c}^{k+q}) end
	(K::Kernal)(k′, k, ψ₁, ψ₂, ψ₃, ψ₄) is used to calculate kernal matrix elements. 
	You can change inputs, but use local operation carefully.
	Do not forget the minus of Kᵈ.

	For the case only refer to q in a specified qgrid, 
	function (K::Kernal)(::Val{:initialize_qgrid}) end
	it's is a supplement to (K::Kernal)(::Val{:initialize}).
	function (K::Kernal)(kq_kindex::Vector{Int}, kΓq_kΓindex::Vector{Int}, q::ReducedCoordinates) end
	for each q, BSEqgrid will not use (K::Kernal)(q::ReducedCoordinates) any more, instead of using this function.
"""
abstract type AbstractKernalInterAction end

function (K::AbstractKernalInterAction)(sym::Symbol, paras...)
	return K(Val(sym), paras...)
end


include("./UJ_trunction/UJ_trunction.jl")
include("./U_trunction/U_trunction.jl")
include("./cluster_UJ_trunction/cluster_UJ_trunction.jl")
