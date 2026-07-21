struct _U_pair_spinaware_Kˣ{
	UT <: DirectInterAction,
	pT <: AbstractReciprocalHoppings,
} <: _spinaware_Kˣ
	nw::Int
	nw_up::Int
	nw_dn::Int
	np::Int
	np_up::Int
	np_dn::Int
	np_ij::Int
	np_up_ij::Int
	np_dn_ij::Int
	nk::Int
	kgrid::RedKgrid
	upindex::Vector{Int}
	dnindex::Vector{Int}
	pair_upindex::Vector{Int}
	pair_dnindex::Vector{Int}
	U::UT
	pair_rh::pT
	Uq::Matrix{ComplexF64}
	Uq_uu::Matrix{ComplexF64}
	Uq_dd::Matrix{ComplexF64}
	Uq_ud::Matrix{ComplexF64}
	Uq_du::Matrix{ComplexF64}
	pairq::Matrix{ComplexF64}
	pairq_uu::Matrix{ComplexF64}
	pairq_dd::Matrix{ComplexF64}
	pairq_ud::Matrix{ComplexF64}
	pairq_du::Matrix{ComplexF64}
	pair::Vector{NTuple{2, Int}}
	pair_ij::Vector{NTuple{2, Int}}
	pair_up::Vector{NTuple{2, Int}}
	pair_up_ij::Vector{NTuple{2, Int}}
	pair_dn::Vector{NTuple{2, Int}}
	pair_dn_ij::Vector{NTuple{2, Int}}
	pair_pR::Vector{ReducedCoordinates{Int}}
	ek_pR::Matrix{ComplexF64}
end
function (V::_U_pair_spinaware_Kˣ)(q::ReducedCoordinates; isΓ::Bool = false)
	if isΓ
		V.U(V.Uq, ReducedCoordinates(0, 0, 0); isΓ = true)
		V.Uq .+= V.U(Val(:taylor_0))
		V.pair_rh(V.pairq, ReducedCoordinates(0, 0, 0))
	else
		V.U(V.Uq, q)
		V.pair_rh(V.pairq, q)
	end
	V.Uq ./= V.nk
	V.Uq_uu .= V.Uq[V.upindex, V.upindex]
	V.Uq_dd .= V.Uq[V.dnindex, V.dnindex]
	V.Uq_ud .= V.Uq[V.upindex, V.dnindex]
	V.Uq_du .= V.Uq[V.dnindex, V.upindex]
	V.pairq ./= V.nk
	V.pairq_uu .= V.pairq[V.pair_upindex, V.pair_upindex]
	V.pairq_dd .= V.pairq[V.pair_dnindex, V.pair_dnindex]
	V.pairq_ud .= V.pairq[V.pair_upindex, V.pair_dnindex]
	V.pairq_du .= V.pairq[V.pair_dnindex, V.pair_upindex]
	return nothing
end
function (V::_U_pair_spinaware_Kˣ)(::Val{:buffer})
	buffer_nw = Vector{ComplexF64}(undef, V.nw)
	buffer_mul_nw = Vector{ComplexF64}(undef, V.nw)
	buffer_np = Vector{ComplexF64}(undef, V.np)
	buffer_mul_np = Vector{ComplexF64}(undef, V.np)
	Up = Vector{ComplexF64}(undef, V.np_ij)
	return buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up
end
function (V::_U_pair_spinaware_Kˣ)(::Val{:buffer_uu})
	buffer_nw = Vector{ComplexF64}(undef, V.nw_up)
	buffer_mul_nw = Vector{ComplexF64}(undef, V.nw_up)
	buffer_np = Vector{ComplexF64}(undef, V.np_up)
	buffer_mul_np = Vector{ComplexF64}(undef, V.np_up)
	Up = Vector{ComplexF64}(undef, V.np_up_ij)
	return buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up
end
function (V::_U_pair_spinaware_Kˣ)(::Val{:buffer_dd})
	buffer_nw = Vector{ComplexF64}(undef, V.nw_dn)
	buffer_mul_nw = Vector{ComplexF64}(undef, V.nw_dn)
	buffer_np = Vector{ComplexF64}(undef, V.np_dn)
	buffer_mul_np = Vector{ComplexF64}(undef, V.np_dn)
	Up = Vector{ComplexF64}(undef, V.np_dn_ij)
	return buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up
end
function (V::_U_pair_spinaware_Kˣ)(::Val{:buffer_uudd})
	buffer_nw_up = Vector{ComplexF64}(undef, V.nw_up)
	buffer_nw_dn = Vector{ComplexF64}(undef, V.nw_dn)
	buffer_mul_nw = Vector{ComplexF64}(undef, V.nw_up)
	buffer_np_up = Vector{ComplexF64}(undef, V.np_up)
	buffer_np_dn = Vector{ComplexF64}(undef, V.np_dn)
	buffer_mul_np = Vector{ComplexF64}(undef, V.np_up)
	Up_up = Vector{ComplexF64}(undef, V.np_up_ij)
	Up_dn = Vector{ComplexF64}(undef, V.np_dn_ij)
	return buffer_nw_up, buffer_nw_dn, buffer_mul_nw, buffer_np_up, buffer_np_dn, buffer_mul_np, Up_up, Up_dn
end
function (V::_U_pair_spinaware_Kˣ)(::Val{:buffer_dduu})
	buffer_nw_up = Vector{ComplexF64}(undef, V.nw_up)
	buffer_nw_dn = Vector{ComplexF64}(undef, V.nw_dn)
	buffer_mul_nw = Vector{ComplexF64}(undef, V.nw_dn)
	buffer_np_up = Vector{ComplexF64}(undef, V.np_up)
	buffer_np_dn = Vector{ComplexF64}(undef, V.np_dn)
	buffer_mul_np = Vector{ComplexF64}(undef, V.np_dn)
	Up_up = Vector{ComplexF64}(undef, V.np_up_ij)
	Up_dn = Vector{ComplexF64}(undef, V.np_dn_ij)
	return buffer_nw_up, buffer_nw_dn, buffer_mul_nw, buffer_np_up, buffer_np_dn, buffer_mul_np, Up_up, Up_dn
end
function (V::_U_pair_spinaware_Kˣ)(ψv′, ψc′, k′, ψv, ψc, k)
	buffer = V(Val(:buffer))
	return V(buffer..., ψv′, ψc′, k′, ψv, ψc, k)
end
function (V::_U_pair_spinaware_Kˣ)(::Val{sym}, ψv′, ψc′, k′, ψv, ψc, k) where {sym}
	buffer = V(Val(Symbol(:buffer_, sym)))
	return V(Val(sym), buffer..., ψv′, ψc′, k′, ψv, ψc, k)
end
function (V::_U_pair_spinaware_Kˣ)(
	buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up,
	ψv′, ψc′, k′, ψv, ψc, k,
)

	pair = V.pair
	pair_ij = V.pair_ij

	buffer_nw .= conj.(ψv)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up[idx] = ψc[i] * buffer_nw[j] # Up
	end
	ek_pR = view(V.ek_pR, :, k)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np[idx] = Up[ij] * ek_pR[R] # Bp
	end

	Bw = view(Up, 1:V.nw)
	mul!(buffer_mul_nw, V.Uq, Bw)
	mul!(buffer_mul_np, V.pairq, buffer_np)

	buffer_nw .= conj.(ψv′)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up[idx] = ψc′[i] * buffer_nw[j] # Up
	end
	ek′_pR = view(V.ek_pR, :, k′)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np[idx] = Up[ij] * ek′_pR[R] # Ap
	end

	Aw = view(Up, 1:V.nw)
	Kˣ = Aw ⋅ buffer_mul_nw
	Kˣ += buffer_np ⋅ buffer_mul_np

	return Kˣ
end
function (V::_U_pair_spinaware_Kˣ)(::Val{:uu},
	buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up,
	ψv′, ψc′, k′, ψv, ψc, k,
)

	pair = V.pair_up
	pair_ij = V.pair_up_ij
	nw = V.nw_up

	buffer_nw .= conj.(ψv)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up[idx] = ψc[i] * buffer_nw[j] # Up
	end
	ek_pR = view(V.ek_pR, :, k)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np[idx] = Up[ij] * ek_pR[R] # Bp
	end

	Bw = view(Up, 1:nw)
	mul!(buffer_mul_nw, V.Uq_uu, Bw)
	mul!(buffer_mul_np, V.pairq_uu, buffer_np)

	buffer_nw .= conj.(ψv′)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up[idx] = ψc′[i] * buffer_nw[j] # Up
	end
	ek′_pR = view(V.ek_pR, :, k′)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np[idx] = Up[ij] * ek′_pR[R] # Ap
	end

	Aw = view(Up, 1:nw)
	Kˣ = Aw ⋅ buffer_mul_nw
	Kˣ += buffer_np ⋅ buffer_mul_np

	return Kˣ
end
function (V::_U_pair_spinaware_Kˣ)(::Val{:dd},
	buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up,
	ψv′, ψc′, k′, ψv, ψc, k,
)

	pair = V.pair_dn
	pair_ij = V.pair_dn_ij
	nw = V.nw_dn

	buffer_nw .= conj.(ψv)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up[idx] = ψc[i] * buffer_nw[j] # Up
	end
	ek_pR = view(V.ek_pR, :, k)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np[idx] = Up[ij] * ek_pR[R] # Bp
	end

	Bw = view(Up, 1:nw)
	mul!(buffer_mul_nw, V.Uq_dd, Bw)
	mul!(buffer_mul_np, V.pairq_dd, buffer_np)

	buffer_nw .= conj.(ψv′)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up[idx] = ψc′[i] * buffer_nw[j] # Up
	end
	ek′_pR = view(V.ek_pR, :, k′)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np[idx] = Up[ij] * ek′_pR[R] # Ap
	end

	Aw = view(Up, 1:nw)
	Kˣ = Aw ⋅ buffer_mul_nw
	Kˣ += buffer_np ⋅ buffer_mul_np

	return Kˣ
end
function (V::_U_pair_spinaware_Kˣ)(::Val{:uudd},
	buffer_nw_up, buffer_nw_dn, buffer_mul_nw, buffer_np_up, buffer_np_dn, buffer_mul_np, Up_up, Up_dn,
	ψv′, ψc′, k′, ψv, ψc, k,
)

	pair = V.pair_dn
	pair_ij = V.pair_dn_ij
	nw = V.nw_dn

	buffer_nw_dn .= conj.(ψv)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up_dn[idx] = ψc[i] * buffer_nw_dn[j] # Up
	end
	ek_pR = view(V.ek_pR, :, k)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np_dn[idx] = Up_dn[ij] * ek_pR[R] # Bp
	end

	Bw = view(Up_dn, 1:nw)
	mul!(buffer_mul_nw, V.Uq_ud, Bw)
	mul!(buffer_mul_np, V.pairq_ud, buffer_np_dn)

	pair = V.pair_up
	pair_ij = V.pair_up_ij
	nw = V.nw_up

	buffer_nw_up .= conj.(ψv′)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up_up[idx] = ψc′[i] * buffer_nw_up[j] # Up
	end
	ek′_pR = view(V.ek_pR, :, k′)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np_up[idx] = Up_up[ij] * ek′_pR[R] # Ap
	end

	Aw = view(Up_up, 1:nw)
	Kˣ = Aw ⋅ buffer_mul_nw
	Kˣ += buffer_np_up ⋅ buffer_mul_np

	return Kˣ
end
function (V::_U_pair_spinaware_Kˣ)(::Val{:dduu},
	buffer_nw_up, buffer_nw_dn, buffer_mul_nw, buffer_np_up, buffer_np_dn, buffer_mul_np, Up_up, Up_dn,
	ψv′, ψc′, k′, ψv, ψc, k,
)

	pair = V.pair_up
	pair_ij = V.pair_up_ij
	nw = V.nw_up

	buffer_nw_up .= conj.(ψv)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up_up[idx] = ψc[i] * buffer_nw_up[j] # Up
	end
	ek_pR = view(V.ek_pR, :, k)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np_up[idx] = Up_up[ij] * ek_pR[R] # Bp
	end

	Bw = view(Up_up, 1:nw)
	mul!(buffer_mul_nw, V.Uq_du, Bw)
	mul!(buffer_mul_np, V.pairq_du, buffer_np_up)

	pair = V.pair_dn
	pair_ij = V.pair_dn_ij
	nw = V.nw_dn

	buffer_nw_dn .= conj.(ψv′)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up_dn[idx] = ψc′[i] * buffer_nw_dn[j] # Up
	end
	ek′_pR = view(V.ek_pR, :, k′)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np_dn[idx] = Up_dn[ij] * ek′_pR[R] # Ap
	end

	Aw = view(Up_dn, 1:nw)
	Kˣ = Aw ⋅ buffer_mul_nw
	Kˣ += buffer_np_dn ⋅ buffer_mul_np

	return Kˣ
end
function _U_pair_spinaware_Kˣ(kgrid::RedKgrid, U::DirectInterAction, pair_int::WanIntijklR; upindex::Vector{Int}, dnindex::Vector{Int})

	nk = length(kgrid)
	nw = numorb(U)
	nw == pair_int.nw || error("Mismatched U and pair interaction!")

	nw_up = length(upindex)
	nw_dn = length(dnindex)

	# pair_up, pair_dn
	pair = pair_int.pair
	pair_upindex = findall(pair) do (iw1, iw2, _)
		iw1 ∈ upindex && iw2 ∈ upindex
	end
	pair_dnindex = findall(pair) do (iw1, iw2, _)
		iw1 ∈ dnindex && iw2 ∈ dnindex
	end

	pair_rh = ReciprocalHoppings(pair_int)

	Uq = Array{ComplexF64}(undef, nw, nw)
	Uq_uu = Matrix{ComplexF64}(undef, nw_up, nw_up)
	Uq_dd = Matrix{ComplexF64}(undef, nw_dn, nw_dn)
	Uq_ud = Matrix{ComplexF64}(undef, nw_up, nw_dn)
	Uq_du = Matrix{ComplexF64}(undef, nw_dn, nw_up)

	np = pair_int.np
	np_up = length(pair_upindex)
	np_dn = length(pair_dnindex)
	pairq = Array{ComplexF64}(undef, np, np)
	pairq_uu = Array{ComplexF64}(undef, np_up, np_up)
	pairq_dd = Array{ComplexF64}(undef, np_dn, np_dn)
	pairq_ud = Array{ComplexF64}(undef, np_up, np_dn)
	pairq_du = Array{ComplexF64}(undef, np_dn, np_up)


	# pair_up, pair_dn. redefine `iw`.
	pair_up = pair[pair_upindex]
	upindex_to_idx = Dict(iw => idx for (idx, iw) in enumerate(upindex))
	pair_up = map(pair_up) do p
		iw1 = upindex_to_idx[p[1]]
		iw2 = upindex_to_idx[p[2]]
		return (iw1, iw2, p[3])
	end
	pair_dn = pair[pair_dnindex]
	dnindex_to_idx = Dict(iw => idx for (idx, iw) in enumerate(dnindex))
	pair_dn = map(pair_dn) do p
		iw1 = dnindex_to_idx[p[1]]
		iw2 = dnindex_to_idx[p[2]]
		return (iw1, iw2, p[3])
	end

	pair_ij = Set{NTuple{2, Int}}()
	for p in pair
		push!(pair_ij, (p[1], p[2]))
	end
	pair_ij = sort(collect(pair_ij))
	pair_ii = [(i, i) for i in 1:nw]
	setdiff!(pair_ij, pair_ii)
	pair_ij = [pair_ii; pair_ij]
	np_ij = length(pair_ij)

	pair_ij_toidx = Dict(p => idx for (idx, p) in enumerate(pair_ij))
	pair = map(pair) do p
		(pair_ij_toidx[(p[1], p[2])], p[3])
	end

	# pair_up
	pair_up_ij = Set{NTuple{2, Int}}()
	for p in pair_up
		push!(pair_up_ij, (p[1], p[2]))
	end
	pair_up_ij = sort(collect(pair_up_ij))
	pair_up_ii = [(i, i) for i in 1:nw_up]
	setdiff!(pair_up_ij, pair_up_ii)
	pair_up_ij = [pair_up_ii; pair_up_ij]
	np_up_ij = length(pair_up_ij)

	pair_up_ij_toidx = Dict(p => idx for (idx, p) in enumerate(pair_up_ij))
	pair_up = map(pair_up) do p
		(pair_up_ij_toidx[(p[1], p[2])], p[3])
	end

	# pair_dn
	pair_dn_ij = Set{NTuple{2, Int}}()
	for p in pair_dn
		push!(pair_dn_ij, (p[1], p[2]))
	end
	pair_dn_ij = sort(collect(pair_dn_ij))
	pair_dn_ii = [(i, i) for i in 1:nw_dn]
	setdiff!(pair_dn_ij, pair_dn_ii)
	pair_dn_ij = [pair_dn_ii; pair_dn_ij]
	np_dn_ij = length(pair_dn_ij)

	pair_dn_ij_toidx = Dict(p => idx for (idx, p) in enumerate(pair_dn_ij))
	pair_dn = map(pair_dn) do p
		(pair_dn_ij_toidx[(p[1], p[2])], p[3])
	end

	ek_pR = [cispi(-2 * (k ⋅ R)) for R in pair_int.pR, k in kgrid]

	return _U_pair_spinaware_Kˣ(nw, nw_up, nw_dn, np, np_up, np_dn, np_ij, np_up_ij, np_dn_ij, nk,
		kgrid, upindex, dnindex, pair_upindex, pair_dnindex, U, pair_rh,
		Uq, Uq_uu, Uq_dd, Uq_ud, Uq_du,
		pairq, pairq_uu, pairq_dd, pairq_ud, pairq_du,
		pair, pair_ij, pair_up, pair_up_ij, pair_dn, pair_dn_ij,
		pair_int.pR, ek_pR)
end
