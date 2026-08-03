struct _U_spinaware_Kˣ{
	UT <: DirectInterAction,
} <: _spinaware_Kˣ
	nw::Int
	nw_up::Int
	nw_dn::Int
	nk::Int
	kgrid::RedKgrid
	upindex::Vector{Int}
	dnindex::Vector{Int}
	U::UT
	Uq::Matrix{ComplexF64}
	Uq_uu::Matrix{ComplexF64}
	Uq_dd::Matrix{ComplexF64}
	Uq_ud::Matrix{ComplexF64}
	Uq_du::Matrix{ComplexF64}
end
function (V::_U_spinaware_Kˣ)(q::ReducedCoordinates; isΓ::Bool = false)
	if isΓ
		V.U(V.Uq, ReducedCoordinates(0, 0, 0); isΓ = true)
		V.Uq .+= V.U(Val(:taylor_0))
	else
		V.U(V.Uq, q)
	end
	V.Uq ./= V.nk
	V.Uq_uu .= V.Uq[V.upindex, V.upindex]
	V.Uq_dd .= V.Uq[V.dnindex, V.dnindex]
	V.Uq_ud .= V.Uq[V.upindex, V.dnindex]
	V.Uq_du .= V.Uq[V.dnindex, V.upindex]
	return nothing
end
function (V::_U_spinaware_Kˣ)(::Val{:buffer})
	buffer_nw = Vector{ComplexF64}(undef, V.nw)
	buffer_mul_nw = Vector{ComplexF64}(undef, V.nw)
	return buffer_nw, buffer_mul_nw
end
function (V::_U_spinaware_Kˣ)(::Val{:buffer_uu})
	buffer_nw = Vector{ComplexF64}(undef, V.nw_up)
	buffer_mul_nw = Vector{ComplexF64}(undef, V.nw_up)
	return buffer_nw, buffer_mul_nw
end
function (V::_U_spinaware_Kˣ)(::Val{:buffer_dd})
	buffer_nw = Vector{ComplexF64}(undef, V.nw_dn)
	buffer_mul_nw = Vector{ComplexF64}(undef, V.nw_dn)
	return buffer_nw, buffer_mul_nw
end
function (V::_U_spinaware_Kˣ)(::Val{:buffer_uudd})
	buffer_nw_up = Vector{ComplexF64}(undef, V.nw_up)
	buffer_nw_dn = Vector{ComplexF64}(undef, V.nw_dn)
	buffer_mul_nw = Vector{ComplexF64}(undef, V.nw_up)
	return buffer_nw_up, buffer_nw_dn, buffer_mul_nw
end
function (V::_U_spinaware_Kˣ)(::Val{:buffer_dduu})
	buffer_nw_up = Vector{ComplexF64}(undef, V.nw_up)
	buffer_nw_dn = Vector{ComplexF64}(undef, V.nw_dn)
	buffer_mul_nw = Vector{ComplexF64}(undef, V.nw_dn)
	return buffer_nw_up, buffer_nw_dn, buffer_mul_nw
end
function (V::_U_spinaware_Kˣ)(ψv′, ψc′, k′, ψv, ψc, k)
	buffer = V(Val(:buffer))
	return V(buffer..., ψv′, ψc′, k′, ψv, ψc, k)
end
function (V::_U_spinaware_Kˣ)(::Val{sym}, ψv′, ψc′, k′, ψv, ψc, k) where {sym}
	buffer = V(Val(Symbol(:buffer_, sym)))
	return V(Val(sym), buffer..., ψv′, ψc′, k′, ψv, ψc, k)
end
function (V::_U_spinaware_Kˣ)(buffer_nw, buffer_mul_nw, ψv′, ψc′, k′, ψv, ψc, k)

	buffer_nw .= ψc .* conj.(ψv) # Bw

	mul!(buffer_mul_nw, V.Uq, buffer_nw)

	buffer_nw .= ψc′ .* conj.(ψv′) # Aw
	Kˣ = buffer_nw ⋅ buffer_mul_nw

	return Kˣ
end
function (V::_U_spinaware_Kˣ)(::Val{:uu}, buffer_nw, buffer_mul_nw, ψv′, ψc′, k′, ψv, ψc, k)

	buffer_nw .= ψc .* conj.(ψv) # Bw

	mul!(buffer_mul_nw, V.Uq_uu, buffer_nw)

	buffer_nw .= ψc′ .* conj.(ψv′) # Aw
	Kˣ = buffer_nw ⋅ buffer_mul_nw

	return Kˣ
end
function (V::_U_spinaware_Kˣ)(::Val{:dd}, buffer_nw, buffer_mul_nw, ψv′, ψc′, k′, ψv, ψc, k)

	buffer_nw .= ψc .* conj.(ψv) # Bw

	mul!(buffer_mul_nw, V.Uq_dd, buffer_nw)

	buffer_nw .= ψc′ .* conj.(ψv′) # Aw
	Kˣ = buffer_nw ⋅ buffer_mul_nw

	return Kˣ
end
function (V::_U_spinaware_Kˣ)(::Val{:uudd}, buffer_nw_up, buffer_nw_dn, buffer_mul_nw, ψv′, ψc′, k′, ψv, ψc, k)

	buffer_nw_dn .= ψc .* conj.(ψv) # Bw

	mul!(buffer_mul_nw, V.Uq_ud, buffer_nw_dn)

	buffer_nw_up .= ψc′ .* conj.(ψv′) # Aw
	Kˣ = buffer_nw_up ⋅ buffer_mul_nw

	return Kˣ
end
function (V::_U_spinaware_Kˣ)(::Val{:dduu}, buffer_nw_up, buffer_nw_dn, buffer_mul_nw, ψv′, ψc′, k′, ψv, ψc, k)

	buffer_nw_up .= ψc .* conj.(ψv) # Bw

	mul!(buffer_mul_nw, V.Uq_du, buffer_nw_up)

	buffer_nw_dn .= ψc′ .* conj.(ψv′) # Aw
	Kˣ = buffer_nw_dn ⋅ buffer_mul_nw

	return Kˣ
end
function _U_spinaware_Kˣ(kgrid::RedKgrid, U::DirectInterAction; upindex::Vector{Int}, dnindex::Vector{Int})

	nw = numorb(U)
	nk = length(kgrid)
	Uq = Array{ComplexF64}(undef, nw, nw)

	nw_up = length(upindex)
	nw_dn = length(dnindex)
	Uq_uu = Matrix{ComplexF64}(undef, nw_up, nw_up)
	Uq_dd = Matrix{ComplexF64}(undef, nw_dn, nw_dn)
	Uq_ud = Matrix{ComplexF64}(undef, nw_up, nw_dn)
	Uq_du = Matrix{ComplexF64}(undef, nw_dn, nw_up)

	return _U_spinaware_Kˣ(nw, nw_up, nw_dn, nk, kgrid, upindex, dnindex, U, Uq, Uq_uu, Uq_dd, Uq_ud, Uq_du)
end
