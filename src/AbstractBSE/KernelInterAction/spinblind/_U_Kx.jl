struct _U_spinblind_Kˣ{
	UT <: DirectInterAction,
} <: _spinblind_Kˣ
	nw::Int
	nk::Int
	kgrid::RedKgrid
	U::UT
	Uq::Matrix{ComplexF64}
end
function (V::_U_spinblind_Kˣ)(q::ReducedCoordinates; isΓ::Bool = false)
	if isΓ
		V.U(V.Uq, ReducedCoordinates(0, 0, 0); isΓ = true)
		V.Uq .+= V.U(Val(:taylor_0))
	else
		V.U(V.Uq, q)
	end
	V.Uq ./= V.nk
	return nothing
end
function (V::_U_spinblind_Kˣ)(::Val{:buffer})
	buffer_nw = Vector{ComplexF64}(undef, V.nw)
	buffer_mul_nw = Vector{ComplexF64}(undef, V.nw)
	return buffer_nw, buffer_mul_nw
end
(V::_U_spinblind_Kˣ)(::Val{:buffer_uu}) = V(Val(:buffer))
(V::_U_spinblind_Kˣ)(::Val{:buffer_dd}) = V(Val(:buffer))
(V::_U_spinblind_Kˣ)(::Val{:buffer_uudd}) = V(Val(:buffer))
(V::_U_spinblind_Kˣ)(::Val{:buffer_dduu}) = V(Val(:buffer))
function (V::_U_spinblind_Kˣ)(ψv′, ψc′, k′, ψv, ψc, k)
	buffer = V(Val(:buffer))
	return V(buffer..., ψv′, ψc′, k′, ψv, ψc, k)

	# Aw = Vector{ComplexF64}(undef, V.nw)
	# Bw = Vector{ComplexF64}(undef, V.nw)
	# buffer_nw = Vector{ComplexF64}(undef, V.nw)

	# Aw .= ψc′ .* conj.(ψv′)
	# Bw .= ψc .* conj.(ψv)

	# mul!(buffer_nw, V.Uq, Bw)
	# Kˣ = Aw ⋅ buffer_nw

	# return Kˣ
end
function (V::_U_spinblind_Kˣ)(buffer_nw, buffer_mul_nw, ψv′, ψc′, k′, ψv, ψc, k)

	buffer_nw .= ψc .* conj.(ψv) # Bw

	mul!(buffer_mul_nw, V.Uq, buffer_nw)

	buffer_nw .= ψc′ .* conj.(ψv′) # Aw
	Kˣ = buffer_nw ⋅ buffer_mul_nw

	return Kˣ
end
function _U_spinblind_Kˣ(kgrid::RedKgrid, U::DirectInterAction)

	nw = numorb(U)
	nk = length(kgrid)
	Uq = Array{ComplexF64}(undef, nw, nw)

	return _U_spinblind_Kˣ(nw, nk, kgrid, U, Uq)
end
