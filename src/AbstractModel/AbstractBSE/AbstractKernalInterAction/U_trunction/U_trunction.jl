
struct _spinless_U_Kᵈ{
	UT <: Number,
	CT <: AbstractLRCorrection,
}
	U::UwithLR{UT, CT}
	Uk::Vector{Matrix{ComplexF64}}
end
struct _spinless_U_Kˣ{
	UT <: Number,
	CT <: AbstractLRCorrection,
}
	U::UwithLR{UT, CT}
	Uq::Matrix{ComplexF64}
end

include("./_spinless_U.jl")
include("./_time_reversal_U.jl")

export Kernal_U

"""
	Kernal_U(kgrid::Union{MonkhorstPack, RedKgrid}, 
		Kᵈ_U::UwithLR, Kˣ_U::UwithLR;
		upindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
		dnindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing
	) -> AbstractKernalInterAction

Creat a object which can be provided to canstruct a BSE model.
- `kgrid`: should be the same as excitonic basis;
- `K_U`: direct Coulomb interaction with long-range correction between electronic wannier basis;
- `up(dn)_index`: the wannier basis index when contains spin.

Its usage is similar to [`Kernal_UJ`](@ref).
!!! note
    For spinful case, we require that the electronic wannier basis satisfy time-reversal symmetry,
	so that we can define which term is `U`.
"""
function Kernal_U(kgrid::MonkhorstPack, paras...; kwards...)
	kgrid = RedKgrid(kgrid)
	return Kernal_U(kgrid, paras...; kwards...)
end
function Kernal_U(kgrid::RedKgrid, Kᵈ_U::UwithLR{UT_d, CT_d}, Kˣ_U::UwithLR{UT_x, CT_x};
	upindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
	dnindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
) where {UT_d, CT_d, UT_x, CT_x}

	nk = length(kgrid)
	norb = Kᵈ_U.norb

	Kᵈ_Uk = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nk]
	W = _spinless_U_Kᵈ(Kᵈ_U, Kᵈ_Uk)

	Kˣ_Uq = Matrix{ComplexF64}(undef, norb, norb)
	V = _spinless_U_Kˣ(Kˣ_U, Kˣ_Uq)


	kgrid_Γ = RedKgrid(MonkhorstPack(kgrid.kgrid_size))
	kgrid_addmap = kgridmap(kgrid_Γ.kdirect, kgrid.kdirect, +)
	kgrid_minusmap = kgridmap(kgrid_Γ.kdirect, kgrid.kdirect, -)

	norb_wispin = 2 * norb
	if isnothing(upindex)
		if isnothing(dnindex)
			return _spinless_U_Kernal(nk, kgrid, kgrid_Γ, kgrid_addmap, kgrid_minusmap, W, V)
		else
			dnindex = collect(Int.(dnindex))
			upindex = setdiff(1:norb_wispin, dnindex)
		end
	else
		upindex = collect(Int.(upindex))
		if isnothing(dnindex)
			dnindex = setdiff(1:norb_wispin, upindex)
		else
			dnindex = collect(Int.(dnindex))
		end
	end
	return _time_reversal_U_Kernal(nk, kgrid, kgrid_Γ, kgrid_addmap, kgrid_minusmap, upindex, dnindex, W, V)
end


function (K::Union{_spinless_U_Kernal, _time_reversal_U_Kernal})(::Val{:initialize})
	Threads.@threads for k in Base.OneTo(K.nk)
		K.W.U(K.W.Uk[k], K.kgrid_Γ[k], K.nk)
		K.W.Uk[k] ./= K.nk
	end
	return nothing
end
function (K::Union{_spinless_U_Kernal, _time_reversal_U_Kernal})(q::ReducedCoordinates)
	#VU(q)
	K.V.U(K.V.Uq, q)
	K.V.Uq ./= K.nk
	return nothing
end
function (K::Union{_spinless_U_Kernal, _time_reversal_U_Kernal})(::Val{:initialize_qgrid})
	return nothing
end
function (K::Union{_spinless_U_Kernal, _time_reversal_U_Kernal})(kq_kindex::Vector{Int}, kΓq_kΓindex::Vector{Int}, q::ReducedCoordinates)
	#VU(q)
	K.V.U(K.V.Uq, q)
	K.V.Uq ./= K.nk
	return nothing
end
