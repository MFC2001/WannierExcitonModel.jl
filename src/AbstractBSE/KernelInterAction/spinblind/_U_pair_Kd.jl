struct _U_pair_spinblind_Kᵈ{
	UT <: DirectInterAction,
	pT <: AbstractReciprocalHoppings,
} <: _spinblind_Kᵈ
	nw::Int
	np::Int
	np_ij::Int
	nk::Int
	kgrid::RedKgrid
	kgrid_Γ::RedKgrid
	kgrid_minusmap::Matrix{Int} #kgrid[k₁] - kgrid[k₂] = kgrid_Γ[kgrid_minusmap[k₁, k₂]]
	U::UT
	pair_rh::pT
	Uk::Array{ComplexF64, 3} #kgrid_Γ
	pairk::Array{ComplexF64, 3} #kgrid_Γ
	pair::Vector{NTuple{2, Int}} #(ij, pR)
	pair_ij::Vector{NTuple{2, Int}} #(i, j)
	pair_pR::Vector{ReducedCoordinates{Int}}
	ek_pR::Matrix{ComplexF64} #[npR, nk]
	eq_pR::Vector{ComplexF64} #[npR]
	ekq_pR::Matrix{ComplexF64} #[npR, nk]
end
function (W::_U_pair_spinblind_Kᵈ)(q::ReducedCoordinates)
	map(enumerate(W.pair_pR)) do (iR, R)
		W.eq_pR[iR] = cispi(-2 * (q ⋅ R))
	end
	W.ekq_pR .= W.ek_pR .* W.eq_pR
	return nothing
end
function (W::_U_pair_spinblind_Kᵈ)(::Val{:buffer})
	buffer_nw = Vector{ComplexF64}(undef, W.nw)
	buffer_mul_nw = Vector{ComplexF64}(undef, W.nw)
	buffer_np = Vector{ComplexF64}(undef, W.np)
	buffer_mul_np = Vector{ComplexF64}(undef, W.np)
	Up = Vector{ComplexF64}(undef, W.np_ij)
	return buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up
end
(W::_U_pair_spinblind_Kᵈ)(::Val{:buffer_uu}) = W(Val(:buffer))
(W::_U_pair_spinblind_Kᵈ)(::Val{:buffer_dd}) = W(Val(:buffer))
(W::_U_pair_spinblind_Kᵈ)(::Val{:buffer_ud}) = W(Val(:buffer))
(W::_U_pair_spinblind_Kᵈ)(::Val{:buffer_du}) = W(Val(:buffer))
function (W::_U_pair_spinblind_Kᵈ)(ψv′, ψc′, k′, ψv, ψc, k)
	buffer = W(Val(:buffer))
	return W(buffer..., ψv′, ψc′, k′, ψv, ψc, k)

	# Aw = Vector{ComplexF64}(undef, W.nw)
	# Bw = Vector{ComplexF64}(undef, W.nw)
	# Ap = Vector{ComplexF64}(undef, W.np)
	# Bp = Vector{ComplexF64}(undef, W.np)
	# buffer_nw = Vector{ComplexF64}(undef, W.nw)
	# buffer_np = Vector{ComplexF64}(undef, W.np)
	# Up = Vector{ComplexF64}(undef, W.np_ij)

	# pair = W.pair
	# pair_ij = W.pair_ij

	# buffer_nw .= conj.(ψc)
	# @inbounds @simd for idx in eachindex(pair_ij)
	# 	(i, j) = pair_ij[idx]
	# 	Up[idx] = ψc′[i] * buffer_nw[j]
	# end
	# Aw .= Up[1:nw]
	# ekq_pR = view(W.ekq_pR, :, k)
	# @inbounds @simd for idx in eachindex(pair)
	# 	(ij, R) = pair[idx]
	# 	Ap[idx] = Up[ij] * ekq_pR[R]
	# end
	# buffer_nw .= conj.(ψv)
	# @inbounds @simd for idx in eachindex(pair_ij)
	# 	(i, j) = pair_ij[idx]
	# 	Up[idx] = ψv′[i] * buffer_nw[j]
	# end
	# Bw .= Up[1:nw]
	# ek_pR = view(W.ek_pR, :, k)
	# @inbounds @simd for idx in eachindex(pair)
	# 	(ij, R) = pair[idx]
	# 	Bp[idx] = Up[ij] * ek_pR[R]
	# end

	# k′k = W.kgrid_minusmap[k′, k]

	# Uk′k = view(W.Uk, :, :, k′k)
	# mul!(buffer_nw, Uk′k, Bw)
	# Kᵈ = Aw ⋅ buffer_nw

	# Pk′k = view(W.pairk, :, :, k′k)
	# mul!(buffer_np, Pk′k, Bp)
	# Kᵈ += Ap ⋅ buffer_np

	# return -Kᵈ
end
function (W::_U_pair_spinblind_Kᵈ)(
	buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up,
	ψv′, ψc′, k′, ψv, ψc, k,
)

	pair = W.pair
	pair_ij = W.pair_ij

	buffer_nw .= conj.(ψv)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up[idx] = ψv′[i] * buffer_nw[j]
	end
	ek_pR = view(W.ek_pR, :, k)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np[idx] = Up[ij] * ek_pR[R] # Bp
	end

	k′k = W.kgrid_minusmap[k′, k]
	Bw = view(Up, 1:W.nw)
	Uk′k = view(W.Uk,:,:,k′k)
	mul!(buffer_mul_nw, Uk′k, Bw)
	Pk′k = view(W.pairk,:,:,k′k)
	mul!(buffer_mul_np, Pk′k, buffer_np)

	buffer_nw .= conj.(ψc)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up[idx] = ψc′[i] * buffer_nw[j]
	end
	ekq_pR = view(W.ekq_pR, :, k)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np[idx] = Up[ij] * ekq_pR[R] # Ap
	end

	Aw = view(Up, 1:W.nw)
	Kᵈ = Aw ⋅ buffer_mul_nw
	Kᵈ += buffer_np ⋅ buffer_mul_np

	return -Kᵈ
end
function _U_pair_spinblind_Kᵈ(kgrid::RedKgrid, U::DirectInterAction, pair_int::WanIntijklR)

	nk = length(kgrid)
	nw = numorb(U)
	nw == pair_int.nw || error("Mismatched U and pair interaction!")

	kgrid_Γ = RedKgrid(MonkhorstPack(kgrid.kgrid_size))
	kgrid_minusmap = kgridmap(kgrid_Γ, kgrid, kgrid, -)
	iΓ = kgrid_minusmap[1, 1]

	pair_rh = ReciprocalHoppings(pair_int)

	Uk = Array{ComplexF64}(undef, nw, nw, nk)
	np = pair_int.np
	pairk = Array{ComplexF64}(undef, np, np, nk)
	Threads.@threads for k in Base.OneTo(nk)
		Uk_view = view(Uk,:,:,k)
		U(Uk_view, kgrid_Γ[k], kgrid.kgrid_size; isΓ = (k == iΓ))
		Uk_view ./= nk
		pairk_view = view(pairk,:,:,k)
		pair_rh(pairk_view, kgrid_Γ[k])
		pairk_view ./= nk
	end

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
	eq_pR = Vector{ComplexF64}(undef, pair_int.npR)
	ekq_pR = similar(ek_pR)

	return _U_pair_spinblind_Kᵈ(nw, np, np_ij, nk, kgrid, kgrid_Γ, kgrid_minusmap, U, pair_rh, Uk, pairk,
		pair, pair_ij, pair_int.pR, ek_pR, eq_pR, ekq_pR)
end
