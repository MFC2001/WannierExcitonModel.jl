struct _U_spinaware_Kᵈ{
	UT <: DirectInterAction,
} <: _spinaware_Kᵈ
	nw::Int
	nw_up::Int
	nw_dn::Int
	nk::Int
	kgrid::RedKgrid
	kgrid_Γ::RedKgrid
	kgrid_minusmap::Matrix{Int} #kgrid[k₁] - kgrid[k₂] = kgrid_Γ[kgrid_minusmap[k₁, k₂]]
	upindex::Vector{Int}
	dnindex::Vector{Int}
	U::UT
	Uk::Array{ComplexF64, 3} #kgrid_Γ
	Uk_uu::Array{ComplexF64, 3}
	Uk_dd::Array{ComplexF64, 3}
	Uk_ud::Array{ComplexF64, 3}
	Uk_du::Array{ComplexF64, 3}
end
(W::_U_spinaware_Kᵈ)(q::ReducedCoordinates) = nothing
function (W::_U_spinaware_Kᵈ)(::Val{:buffer})
	buffer_nw = Vector{ComplexF64}(undef, W.nw)
	buffer_mul_nw = Vector{ComplexF64}(undef, W.nw)
	return buffer_nw, buffer_mul_nw
end
function (W::_U_spinaware_Kᵈ)(::Val{:buffer_uu})
	buffer_nw = Vector{ComplexF64}(undef, W.nw_up)
	buffer_mul_nw = Vector{ComplexF64}(undef, W.nw_up)
	return buffer_nw, buffer_mul_nw
end
function (W::_U_spinaware_Kᵈ)(::Val{:buffer_dd})
	buffer_nw = Vector{ComplexF64}(undef, W.nw_dn)
	buffer_mul_nw = Vector{ComplexF64}(undef, W.nw_dn)
	return buffer_nw, buffer_mul_nw
end
function (W::_U_spinaware_Kᵈ)(::Val{:buffer_ud})
	buffer_nw_up = Vector{ComplexF64}(undef, W.nw_up)
	buffer_nw_dn = Vector{ComplexF64}(undef, W.nw_dn)
	buffer_mul_nw = Vector{ComplexF64}(undef, W.nw_dn)
	return buffer_nw_up, buffer_nw_dn, buffer_mul_nw
end
function (W::_U_spinaware_Kᵈ)(::Val{:buffer_du})
	buffer_nw_up = Vector{ComplexF64}(undef, W.nw_up)
	buffer_nw_dn = Vector{ComplexF64}(undef, W.nw_dn)
	buffer_mul_nw = Vector{ComplexF64}(undef, W.nw_up)
	return buffer_nw_up, buffer_nw_dn, buffer_mul_nw
end
function (W::_U_spinaware_Kᵈ)(ψv′, ψc′, k′, ψv, ψc, k)
	buffer = W(Val(:buffer))
	return W(buffer..., ψv′, ψc′, k′, ψv, ψc, k)
end
function (W::_U_spinaware_Kᵈ)(::Val{sym}, ψv′, ψc′, k′, ψv, ψc, k) where {sym}
	buffer = W(Val(Symbol(:buffer_, sym)))
	return W(Val(sym), buffer..., ψv′, ψc′, k′, ψv, ψc, k)
end
function (W::_U_spinaware_Kᵈ)(buffer_nw, buffer_mul_nw, ψv′, ψc′, k′, ψv, ψc, k)

	buffer_nw .= ψv′ .* conj.(ψv) # Bw

	Uk′k = view(W.Uk,:,:,W.kgrid_minusmap[k′, k])
	mul!(buffer_mul_nw, Uk′k, buffer_nw)

	buffer_nw .= ψc′ .* conj.(ψc) # Aw
	Kᵈ = buffer_nw ⋅ buffer_mul_nw

	return -Kᵈ
end
function (W::_U_spinaware_Kᵈ)(::Val{:uu}, buffer_nw, buffer_mul_nw, ψv′, ψc′, k′, ψv, ψc, k)

	buffer_nw .= ψv′ .* conj.(ψv) # Bw

	Uk′k = view(W.Uk_uu,:,:,W.kgrid_minusmap[k′, k])
	mul!(buffer_mul_nw, Uk′k, buffer_nw)

	buffer_nw .= ψc′ .* conj.(ψc) # Aw
	Kᵈ = buffer_nw ⋅ buffer_mul_nw

	return -Kᵈ
end
function (W::_U_spinaware_Kᵈ)(::Val{:dd}, buffer_nw, buffer_mul_nw, ψv′, ψc′, k′, ψv, ψc, k)

	buffer_nw .= ψv′ .* conj.(ψv) # Bw

	Uk′k = view(W.Uk_dd,:,:,W.kgrid_minusmap[k′, k])
	mul!(buffer_mul_nw, Uk′k, buffer_nw)

	buffer_nw .= ψc′ .* conj.(ψc) # Aw
	Kᵈ = buffer_nw ⋅ buffer_mul_nw

	return -Kᵈ
end
function (W::_U_spinaware_Kᵈ)(::Val{:ud}, buffer_nw_up, buffer_nw_dn, buffer_mul_nw, ψv′, ψc′, k′, ψv, ψc, k)

	buffer_nw_up .= ψv′ .* conj.(ψv) # Bw

	Uk′k = view(W.Uk_du,:,:,W.kgrid_minusmap[k′, k])
	mul!(buffer_mul_nw, Uk′k, buffer_nw_up)

	buffer_nw_dn .= ψc′ .* conj.(ψc) # Aw
	Kᵈ = buffer_nw_dn ⋅ buffer_mul_nw

	return -Kᵈ
end
function (W::_U_spinaware_Kᵈ)(::Val{:du}, buffer_nw_up, buffer_nw_dn, buffer_mul_nw, ψv′, ψc′, k′, ψv, ψc, k)

	buffer_nw_dn .= ψv′ .* conj.(ψv) # Bw

	Uk′k = view(W.Uk_ud,:,:,W.kgrid_minusmap[k′, k])
	mul!(buffer_mul_nw, Uk′k, buffer_nw_dn)

	buffer_nw_up .= ψc′ .* conj.(ψc) # Aw
	Kᵈ = buffer_nw_up ⋅ buffer_mul_nw

	return -Kᵈ
end
function _U_spinaware_Kᵈ(kgrid::RedKgrid, U::DirectInterAction; upindex::Vector{Int}, dnindex::Vector{Int})

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

	Uk_uu = Uk[upindex, upindex, :]
	Uk_dd = Uk[dnindex, dnindex, :]
	Uk_ud = Uk[upindex, dnindex, :]
	Uk_du = Uk[dnindex, upindex, :]

	nw_up = length(upindex)
	nw_dn = length(dnindex)

	return _U_spinaware_Kᵈ(nw, nw_up, nw_dn, nk, kgrid, kgrid_Γ, kgrid_minusmap, upindex, dnindex,
		U, Uk, Uk_uu, Uk_dd, Uk_ud, Uk_du)
end
