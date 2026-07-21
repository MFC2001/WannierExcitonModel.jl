struct _U_J_spinblind_Kᵈ{
	UT <: DirectInterAction,
	J¹T <: AbstractReciprocalHoppings,
	J²T <: AbstractReciprocalHoppings,
} <: _spinblind_Kᵈ
	nw::Int
	nk::Int
	kgrid::RedKgrid
	kgrid_Γ::RedKgrid
	kgrid_minusmap::Matrix{Int} #kgrid[k₁] - kgrid[k₂] = kgrid_Γ[kgrid_minusmap[k₁, k₂]]
	kgrid_addmap::Matrix{Int} #kgrid[k₁] + kgrid[k₂] = kgrid_Γ[kgrid_addmap[k₁, k₂]]
	U::UT
	J¹::J¹T
	J²::J²T
	Uk::Array{ComplexF64, 3} # kgrid_Γ
	J¹q::Matrix{ComplexF64}
	J²kq::Array{ComplexF64, 3} # kgrid_Γ + q
end
function (W::_U_J_spinblind_Kᵈ)(q::ReducedCoordinates)
	W.J¹(W.J¹q, q)
	W.J¹q ./= W.nk
	Threads.@threads for k in Base.OneTo(W.nk)
		kq = W.kgrid_Γ[k] + q
		J²kq = view(W.J²kq,:,:,k)
		W.J²(J²kq, kq)
		J²kq ./= W.nk
	end
	return nothing
end
function (W::_U_J_spinblind_Kᵈ)(::Val{:buffer})
	buffer_nw = Vector{ComplexF64}(undef, W.nw)
	buffer_mul_nw = Vector{ComplexF64}(undef, W.nw)
	return buffer_nw, buffer_mul_nw
end
(W::_U_J_spinblind_Kᵈ)(::Val{:buffer_uu}) = W(Val(:buffer))
(W::_U_J_spinblind_Kᵈ)(::Val{:buffer_dd}) = W(Val(:buffer))
(W::_U_J_spinblind_Kᵈ)(::Val{:buffer_ud}) = W(Val(:buffer))
(W::_U_J_spinblind_Kᵈ)(::Val{:buffer_du}) = W(Val(:buffer))
function (W::_U_J_spinblind_Kᵈ)(ψv′, ψc′, k′, ψv, ψc, k)
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

	# Aw .= ψc′ .* conj.(ψv′)
	# Bw .= ψc .* conj.(ψv)

	# mul!(buffer_nw, W.J¹q, Bw)
	# Kᵈ += Aw ⋅ buffer_nw

	# Aw .= ψc′ .* ψv
	# Bw .= ψc .* ψv′

	# J²k′kq = view(W.J²kq, :, :, W.kgrid_addmap[k′, k])
	# mul!(buffer_nw, J²k′kq, Bw)
	# Kᵈ += Aw ⋅ buffer_nw

	# return -Kᵈ
end
function (W::_U_J_spinblind_Kᵈ)(buffer_nw, buffer_mul_nw, ψv′, ψc′, k′, ψv, ψc, k)

	# Uk′k
	buffer_nw .= ψv′ .* conj.(ψv) # Bw
	Uk′k = view(W.Uk,:,:,W.kgrid_minusmap[k′, k])
	mul!(buffer_mul_nw, Uk′k, buffer_nw)

	buffer_nw .= ψc′ .* conj.(ψc) # Aw
	Kᵈ = buffer_nw ⋅ buffer_mul_nw

	# J¹q
	buffer_nw .= ψc .* conj.(ψv) # Bw
	mul!(buffer_mul_nw, W.J¹q, buffer_nw)

	buffer_nw .= ψc′ .* conj.(ψv′) # Aw
	Kᵈ += buffer_nw ⋅ buffer_mul_nw

	# J²k′kq
	buffer_nw .= ψc .* ψv′ # Bw
	J²k′kq = view(W.J²kq,:,:,W.kgrid_addmap[k′, k])
	mul!(buffer_mul_nw, J²k′kq, buffer_nw)

	buffer_nw .= ψc′ .* ψv # Aw
	Kᵈ += buffer_nw ⋅ buffer_mul_nw

	return -Kᵈ
end
function _U_J_spinblind_Kᵈ(kgrid::RedKgrid, U::DirectInterAction, J¹::AbstractReciprocalHoppings, J²::AbstractReciprocalHoppings)

	nk = length(kgrid)
	nw = numorb(U)
	(nw == numorb(J¹) && nw == numorb(J²)) || error("Mismatched U and pair interaction!")

	kgrid_Γ = RedKgrid(MonkhorstPack(kgrid.kgrid_size))
	kgrid_minusmap = kgridmap(kgrid_Γ, kgrid, kgrid, -)
	kgrid_addmap = kgridmap(kgrid_Γ, kgrid, kgrid, +)
	iΓ = kgrid_minusmap[1, 1]

	Uk = Array{ComplexF64}(undef, nw, nw, nk)
	Threads.@threads for k in Base.OneTo(nk)
		Uk_view = view(Uk,:,:,k)
		U(Uk_view, kgrid_Γ[k], kgrid.kgrid_size; isΓ = (k == iΓ))
		Uk_view ./= nk
	end

	J¹q = Matrix{ComplexF64}(undef, nw, nw)
	J²kq = Array{ComplexF64}(undef, nw, nw, nk)

	return _U_J_spinblind_Kᵈ(nw, nk, kgrid, kgrid_Γ, kgrid_minusmap, kgrid_addmap,
		U, J¹, J², Uk, J¹q, J²kq)
end
