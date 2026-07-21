struct _U_J_spinblind_Kˣ{
	UT <: DirectInterAction,
	J¹T <: AbstractReciprocalHoppings,
	J²T <: AbstractReciprocalHoppings,
} <: _spinblind_Kˣ
	nw::Int
	nk::Int
	kgrid::RedKgrid
	kgrid_Γ::RedKgrid
	kgrid_minusmap::Matrix{Int} #kgrid[k₁] - kgrid[k₂] = kgrid_Γ[kgrid_minusmap[k₁, k₂]]
	kgrid_addmap::Matrix{Int} #kgrid[k₁] + kgrid[k₂] = kgrid_Γ[kgrid_addmap[k₁, k₂]]
	U::UT
	J¹::J¹T
	J²::J²T
	Uq::Matrix{ComplexF64}
	J¹k::Array{ComplexF64, 3}
	J²kq::Array{ComplexF64, 3}
end
function (V::_U_J_spinblind_Kˣ)(q::ReducedCoordinates; isΓ::Bool = false)
	if isΓ
		V.U(V.Uq, ReducedCoordinates(0, 0, 0); isΓ = true)
		V.Uq .+= V.U(Val(:taylor_0))
		Threads.@threads for k in Base.OneTo(V.nk)
			kq = V.kgrid_Γ[k]
			J²kq = view(V.J²kq,:,:,k)
			V.J²(J²kq, kq)
			J²kq ./= V.nk
		end
	else
		V.U(V.Uq, q)
		Threads.@threads for k in Base.OneTo(V.nk)
			kq = V.kgrid_Γ[k] + q
			J²kq = view(V.J²kq,:,:,k)
			V.J²(J²kq, kq)
			J²kq ./= V.nk
		end
	end
	V.Uq ./= V.nk
	return nothing
end
function (V::_U_J_spinblind_Kˣ)(::Val{:buffer})
	buffer_nw = Vector{ComplexF64}(undef, V.nw)
	buffer_mul_nw = Vector{ComplexF64}(undef, V.nw)
	return buffer_nw, buffer_mul_nw
end
(V::_U_J_spinblind_Kˣ)(::Val{:buffer_uu}) = V(Val(:buffer))
(V::_U_J_spinblind_Kˣ)(::Val{:buffer_dd}) = V(Val(:buffer))
(V::_U_J_spinblind_Kˣ)(::Val{:buffer_uudd}) = V(Val(:buffer))
(V::_U_J_spinblind_Kˣ)(::Val{:buffer_dduu}) = V(Val(:buffer))
function (V::_U_J_spinblind_Kˣ)(ψv′, ψc′, k′, ψv, ψc, k)
	buffer = V(Val(:buffer))
	return V(buffer..., ψv′, ψc′, k′, ψv, ψc, k)

	# Aw = Vector{ComplexF64}(undef, V.nw)
	# Bw = Vector{ComplexF64}(undef, V.nw)
	# buffer_nw = Vector{ComplexF64}(undef, V.nw)

	# Aw .= ψc′ .* conj.(ψv′)
	# Bw .= ψc .* conj.(ψv)

	# mul!(buffer_nw, V.Uq, Bw)
	# Kˣ = Aw ⋅ buffer_nw

	# Aw .= ψc′ .* conj.(ψc)
	# Bw .= ψv′ .* conj.(ψv)

	# J¹k′k = view(V.J¹k, :, :, V.kgrid_minusmap[k′, k])
	# mul!(buffer_nw, J¹k′k, Bw)
	# Kˣ += Aw ⋅ buffer_nw

	# Aw .= ψc′ .* ψv
	# Bw .= ψc .* ψv′

	# J²k′kq = view(V.J²kq, :, :, V.kgrid_addmap[k′, k])
	# mul!(buffer_nw, J²k′kq, Bw)
	# Kˣ += Aw ⋅ buffer_nw

	# return Kˣ
end
function (V::_U_J_spinblind_Kˣ)(buffer_nw, buffer_mul_nw, ψv′, ψc′, k′, ψv, ψc, k)

	# Uq
	buffer_nw .= ψc .* conj.(ψv) # Bw
	mul!(buffer_mul_nw, V.Uq, buffer_nw)

	buffer_nw .= ψc′ .* conj.(ψv′) # Aw
	Kˣ = buffer_nw ⋅ buffer_mul_nw

	# J¹k′k
	buffer_nw .= ψv′ .* conj.(ψv) # Bw
	J¹k′k = view(V.J¹k,:,:,V.kgrid_minusmap[k′, k])
	mul!(buffer_mul_nw, J¹k′k, buffer_nw)

	buffer_nw .= ψc′ .* conj.(ψc) # Aw
	Kˣ += buffer_nw ⋅ buffer_mul_nw

	# J²k′kq
	buffer_nw .= ψc .* ψv′ # Bw
	J²k′kq = view(V.J²kq,:,:,V.kgrid_addmap[k′, k])
	mul!(buffer_mul_nw, J²k′kq, buffer_nw)

	buffer_nw .= ψc′ .* ψv # Aw
	Kˣ += buffer_nw ⋅ buffer_mul_nw

	return Kˣ
end
function _U_J_spinblind_Kˣ(kgrid::RedKgrid, U::DirectInterAction, J¹::AbstractReciprocalHoppings, J²::AbstractReciprocalHoppings)

	nk = length(kgrid)
	nw = numorb(U)
	(nw == numorb(J¹) && nw == numorb(J²)) || error("Mismatched U and pair interaction!")

	kgrid_Γ = RedKgrid(MonkhorstPack(kgrid.kgrid_size))
	kgrid_minusmap = kgridmap(kgrid_Γ, kgrid, kgrid, -)
	kgrid_addmap = kgridmap(kgrid_Γ, kgrid, kgrid, +)

	J¹k = Array{ComplexF64}(undef, nw, nw, nk)
	Threads.@threads for k in Base.OneTo(nk)
		J¹k_view = view(J¹k,:,:,k)
		J¹(J¹k_view, kgrid_Γ[k])
		J¹k_view ./= nk
	end

	Uq = Matrix{ComplexF64}(undef, nw, nw)
	J²kq = Array{ComplexF64}(undef, nw, nw, nk)

	return _U_J_spinblind_Kˣ(nw, nk, kgrid, kgrid_Γ, kgrid_minusmap, kgrid_addmap,
		U, J¹, J², Uq, J¹k, J²kq)
end
