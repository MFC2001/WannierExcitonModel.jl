export BSE

abstract type AbstractBSE <: AbstractModel end

# TODO: 应效仿WannierInterAction.jl预先计算密度矩阵的方式，预先计算pair中涉及到的密度矩阵元(此计算仅涉及到简单的乘法和加法运算)
# TODO: 平面波基组下密度矩阵元维度为nG，这里为np；这里np越大，内存开销越大，但同时节省的计算量也越大。
# TODO: 这里通过一点点基组顺序的设计，以及同时只存储需要的密度矩阵元，可以将内存开销降低至极小，约为1e2MB量级。
# TODO: 主要是Kd需要一些设计，Kx的计算很简单直接。

# TODO: 调整kernel的计算，目前的计算瓶颈应该主要是kernel的计算中关于pair的部分，目前的方式多次读取pair相互作用矩阵会导致内存带宽受限；
# TODO: 矩阵维度超过700左右时，就会明显地产生内存带宽的问题。
# TODO: with or without TDA，写两套，分析Kd和Kx，尝试兼容。

include("./basis/basis.jl")
include("./_BSE_Gammadata.jl")
include("./wiTDA/wiTDA.jl")
include("./woTDA/woTDA.jl")

numatom(bse::AbstractBSE) = length(bse.atom_location)
WannierInterActionBase.numorb(bse::AbstractBSE) = length(bse.orb_location)

(bse::AbstractBSE)(sym::Symbol, args...; kwargs...) = bse(Val(sym), args...; kwargs...)
(bse::AbstractBSE)(q::AbstractVector{<:Real}; kwargs...) = bse(ReducedCoordinates(q); kwargs...)

"""
	BSE(TB::AbstractTightBindModel, Kernel::KernelInterAction, sym::Symbol;
		kgrid::Union{MonkhorstPack, RedKgrid} = MonkhorstPack([1, 1, 1]), 
		v::AbstractVector{<:Integer}, c::AbstractVector{<:Integer}, 
		scissor::Real = 0, 
		isqgrid::Bool = false) -> AbstractBSE

Create a BSE model, the type depends on `sym`.
- `TB::AbstractTightBindModel`: contains crystal structure and electronic band structure;
- `Kernel::KernelInterAction`: describes how to calculate `Kᵈ` and `Kˣ`;
- `sym::Symbol`: can be set as `:SU2`, `:general`, `:cluster_SU2`, `:cluster_general`;
- `kgrid::Union{MonkhorstPack, RedKgrid}`: defines the electronic kgrid used to construct the basis of exciton;
- `v::AbstractVector{<:Integer}`: the valence band index used to construct the basis of exciton;
- `c::AbstractVector{<:Integer}`: the conduction band index used to construct the basis of exciton;
- `scissor::Real`: scissor operator used to provide electronic band gap correction;
- `isqgrid::Bool`: when all the excitonic momentum needed can be regarded as the momentum difference in the kgrid.
"""
function BSE(sym::Symbol, args...; TDA = true, kwargs...)
	return BSE(Val(sym), Val(TDA), args...; kwargs...)
end

#From Rohlfing and Louie, we need to control the arbitrary phase of electron wavefunction.
# Acctually, the random phase of electron wavefunction does not change the excitonic result;
# But a definite phase is helpful to get the exact electron wavefunction anytime.

"""
这里对Γ的预处理，主要是为了明确Kˣ的计算方式，其由isΓ决定.
如果isΓ, 则q仅代表方向。
"""
function _BSE_preprocess_Γ(q, isΓ, period)
	if iszero(q) || all(isinteger, q)
		# q等价于Γ.
		np = count(period)
		if np ≠ 3
			# only 3D case, we need to add a term which depends on the direction of q.
			return ReducedCoordinates(0, 0, 0), true
		else
			# 3D case
			if isΓ
				if iszero(q)
					@info "A zero vector `q` is unacceptable for Γ calculation in 3D system."
				else
					# 指定isΓ，此时q代表了方向。
					return ReducedCoordinates(q), true
				end
			else
				@info "A zero vector or an integer vector `q` is unacceptable for non-Γ calculation in 3D system, set `isΓ` as true."
			end
			@info "set ̂q as [0, 0, 1]."
			return ReducedCoordinates(0, 0, 1), true
		end
	else
		# q不等价于Γ.
		return q, isΓ
	end
end
