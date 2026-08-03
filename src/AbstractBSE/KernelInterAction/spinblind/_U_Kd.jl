struct _U_spinblind_Kᵈ{
	UT <: DirectInterAction,
} <: _spinblind_Kᵈ
	nw::Int
	nk::Int
	kgrid::RedKgrid
	kgrid_Γ::RedKgrid
	kgrid_minusmap::Matrix{Int} #kgrid[k₁] - kgrid[k₂] = kgrid_Γ[kgrid_minusmap[k₁, k₂]]
	U::UT
	Uk::Array{ComplexF64, 3} #kgrid_Γ
end
(W::_U_spinblind_Kᵈ)(q::ReducedCoordinates) = nothing
function (W::_U_spinblind_Kᵈ)(::Val{:buffer})
	buffer_nw = Vector{ComplexF64}(undef, W.nw)
	buffer_mul_nw = Vector{ComplexF64}(undef, W.nw)
	return buffer_nw, buffer_mul_nw
end
(W::_U_spinblind_Kᵈ)(::Val{:buffer_uu}) = W(Val(:buffer))
(W::_U_spinblind_Kᵈ)(::Val{:buffer_dd}) = W(Val(:buffer))
(W::_U_spinblind_Kᵈ)(::Val{:buffer_ud}) = W(Val(:buffer))
(W::_U_spinblind_Kᵈ)(::Val{:buffer_du}) = W(Val(:buffer))
function (W::_U_spinblind_Kᵈ)(ψv′, ψc′, k′, ψv, ψc, k)
	buffer = W(Val(:buffer))
	return W(buffer..., ψv′, ψc′, k′, ψv, ψc, k)

	# Aw = Vector{ComplexF64}(undef, W.nw)
	# Bw = Vector{ComplexF64}(undef, W.nw)
	# buffer_nw = Vector{ComplexF64}(undef, W.nw)

	# Aw .= ψc′ .* conj.(ψc)
	# Bw .= ψv′ .* conj.(ψv)

	# Uk′k = view(W.Uk, :, :, W.kgrid_minusmap[k′, k])
	# mul!(buffer_nw, Uk′k, Bw)
	# Kᵈ = Aw ⋅ buffer_nw

	# return -Kᵈ
end
function (W::_U_spinblind_Kᵈ)(buffer_nw, buffer_mul_nw, ψv′, ψc′, k′, ψv, ψc, k)

	buffer_nw .= ψv′ .* conj.(ψv) # Bw

	Uk′k = view(W.Uk,:,:,W.kgrid_minusmap[k′, k])
	mul!(buffer_mul_nw, Uk′k, buffer_nw)

	buffer_nw .= ψc′ .* conj.(ψc) # Aw
	Kᵈ = buffer_nw ⋅ buffer_mul_nw

	return -Kᵈ
end
function _U_spinblind_Kᵈ(kgrid::RedKgrid, U::DirectInterAction)

	kgrid_Γ = RedKgrid(MonkhorstPack(kgrid.kgrid_size))
	kgrid_minusmap = kgridmap(kgrid_Γ, kgrid, kgrid, -)
	iΓ = kgrid_minusmap[1, 1]

	nw = numorb(U)
	nk = length(kgrid)
	Uk = Array{ComplexF64}(undef, nw, nw, nk)
	Threads.@threads for k in Base.OneTo(nk)
		Uk_view = view(Uk,:,:,k)
		U(Uk_view, kgrid_Γ[k], kgrid.kgrid_size; isΓ = (k == iΓ))
		Uk_view ./= nk
	end

	return _U_spinblind_Kᵈ(nw, nk, kgrid, kgrid_Γ, kgrid_minusmap, U, Uk)
end
