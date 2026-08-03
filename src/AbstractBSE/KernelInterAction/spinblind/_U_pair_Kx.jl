struct _U_pair_spinblind_Kˣ{
	UT <: DirectInterAction,
	pT <: AbstractReciprocalHoppings,
} <: _spinblind_Kˣ
	nw::Int
	np::Int
	np_ij::Int
	nk::Int
	kgrid::RedKgrid
	U::UT
	pair_rh::pT
	Uq::Matrix{ComplexF64}
	pairq::Matrix{ComplexF64}
	pair::Vector{NTuple{2, Int}}
	pair_ij::Vector{NTuple{2, Int}}
	pair_pR::Vector{ReducedCoordinates{Int}}
	ek_pR::Matrix{ComplexF64}
end
function (V::_U_pair_spinblind_Kˣ)(q::ReducedCoordinates; isΓ::Bool = false)
	if isΓ
		V.U(V.Uq, ReducedCoordinates(0, 0, 0); isΓ = true)
		V.Uq .+= V.U(Val(:taylor_0))
		V.pair_rh(V.pairq, ReducedCoordinates(0, 0, 0))
	else
		V.U(V.Uq, q)
		V.pair_rh(V.pairq, q)
	end
	V.Uq ./= V.nk
	V.pairq ./= V.nk
	return nothing
end
function (V::_U_pair_spinblind_Kˣ)(::Val{:buffer})
	buffer_nw = Vector{ComplexF64}(undef, V.nw)
	buffer_mul_nw = Vector{ComplexF64}(undef, V.nw)
	buffer_np = Vector{ComplexF64}(undef, V.np)
	buffer_mul_np = Vector{ComplexF64}(undef, V.np)
	Up = Vector{ComplexF64}(undef, V.np_ij)
	return buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up
end
(V::_U_pair_spinblind_Kˣ)(::Val{:buffer_uu}) = V(Val(:buffer))
(V::_U_pair_spinblind_Kˣ)(::Val{:buffer_dd}) = V(Val(:buffer))
(V::_U_pair_spinblind_Kˣ)(::Val{:buffer_uudd}) = V(Val(:buffer))
(V::_U_pair_spinblind_Kˣ)(::Val{:buffer_dduu}) = V(Val(:buffer))
function (V::_U_pair_spinblind_Kˣ)(ψv′, ψc′, k′, ψv, ψc, k)

	buffer = V(Val(:buffer))
	return V(buffer..., ψv′, ψc′, k′, ψv, ψc, k)

	# Aw = Vector{ComplexF64}(undef, V.nw)
	# Bw = Vector{ComplexF64}(undef, V.nw)
	# Ap = Vector{ComplexF64}(undef, V.np)
	# Bp = Vector{ComplexF64}(undef, V.np)
	# buffer_nw = Vector{ComplexF64}(undef, V.nw)
	# buffer_np = Vector{ComplexF64}(undef, V.np)
	# Up = Vector{ComplexF64}(undef, V.np_ij)

	# pair = V.pair
	# pair_ij = V.pair_ij

	# buffer_nw .= conj.(ψv′)
	# @inbounds @simd for idx in eachindex(pair_ij)
	# 	(i, j) = pair_ij[idx]
	# 	Up[idx] = ψc′[i] * buffer_nw[j]
	# end
	# Aw .= Up[1:nw]
	# ek′_pR = view(V.ek_pR, :, k′)
	# @inbounds @simd for idx in eachindex(pair)
	# 	(ij, R) = pair[idx]
	# 	Ap[idx] = Up[ij] * ek′_pR[R]
	# end
	# buffer_nw .= conj.(ψv)
	# @inbounds @simd for idx in eachindex(pair_ij)
	# 	(i, j) = pair_ij[idx]
	# 	Up[idx] = ψc[i] * buffer_nw[j]
	# end
	# Bw .= Up[1:nw]
	# ek_pR = view(V.ek_pR, :, k)
	# @inbounds @simd for idx in eachindex(pair)
	# 	(ij, R) = pair[idx]
	# 	Bp[idx] = Up[ij] * ek_pR[R]
	# end

	# mul!(buffer_nw, V.Uq, Bw)
	# Kˣ = Aw ⋅ buffer_nw

	# mul!(buffer_np, V.pairq, Bp)
	# Kˣ += Ap ⋅ buffer_np

	# return Kˣ
end
function (V::_U_pair_spinblind_Kˣ)(
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
function _U_pair_spinblind_Kˣ(kgrid::RedKgrid, U::DirectInterAction, pair_int::WanIntijklR)

	nk = length(kgrid)
	nw = numorb(U)
	nw == pair_int.nw || error("Mismatched U and pair interaction!")

	pair_rh = ReciprocalHoppings(pair_int)

	Uq = Array{ComplexF64}(undef, nw, nw)
	np = pair_int.np
	pairq = Array{ComplexF64}(undef, np, np)

	pair = pair_int.pair
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

	ek_pR = [cispi(-2 * (k ⋅ R)) for R in pair_int.pR, k in kgrid]

	return _U_pair_spinblind_Kˣ(nw, np, np_ij, nk, kgrid, U, pair_rh, Uq, pairq,
		pair, pair_ij, pair_int.pR, ek_pR)
end
