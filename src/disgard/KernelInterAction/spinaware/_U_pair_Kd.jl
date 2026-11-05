struct _U_pair_spinaware_Kᵈ{
	UT <: DirectInterAction,
	pT <: AbstractReciprocalHoppings,
} <: _spinaware_Kᵈ
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
	kgrid_Γ::RedKgrid
	kgrid_minusmap::Matrix{Int} #kgrid[k₁] - kgrid[k₂] = kgrid_Γ[kgrid_minusmap[k₁, k₂]]
	upindex::Vector{Int}
	dnindex::Vector{Int}
	pair_upindex::Vector{Int}
	pair_dnindex::Vector{Int}
	U::UT
	pair_rh::pT
	Uk::Array{ComplexF64, 3} #kgrid_Γ
	Uk_uu::Array{ComplexF64, 3}
	Uk_dd::Array{ComplexF64, 3}
	Uk_ud::Array{ComplexF64, 3}
	Uk_du::Array{ComplexF64, 3}
	pairk::Array{ComplexF64, 3} #kgrid_Γ
	pairk_uu::Array{ComplexF64, 3}
	pairk_dd::Array{ComplexF64, 3}
	pairk_ud::Array{ComplexF64, 3}
	pairk_du::Array{ComplexF64, 3}
	pair::Vector{NTuple{2, Int}} #(ij, pR)
	pair_ij::Vector{NTuple{2, Int}} #(i, j)
	pair_up::Vector{NTuple{2, Int}}
	pair_up_ij::Vector{NTuple{2, Int}}
	pair_dn::Vector{NTuple{2, Int}}
	pair_dn_ij::Vector{NTuple{2, Int}}
	pair_pR::Vector{ReducedCoordinates{Int}}
	ek_pR::Matrix{ComplexF64} #[npR, nk]
	eq_pR::Vector{ComplexF64} #[npR]
	ekq_pR::Matrix{ComplexF64} #[npR, nk]
end
function (W::_U_pair_spinaware_Kᵈ)(q::ReducedCoordinates)
	map(enumerate(W.pair_pR)) do (iR, R)
		W.eq_pR[iR] = cispi(-2 * (q ⋅ R))
	end
	W.ekq_pR .= W.ek_pR .* W.eq_pR
	return nothing
end
function (W::_U_pair_spinaware_Kᵈ)(::Val{:buffer})
	buffer_nw = Vector{ComplexF64}(undef, W.nw)
	buffer_mul_nw = Vector{ComplexF64}(undef, W.nw)
	buffer_np = Vector{ComplexF64}(undef, W.np)
	buffer_mul_np = Vector{ComplexF64}(undef, W.np)
	Up = Vector{ComplexF64}(undef, W.np_ij)
	return buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up
end
function (W::_U_pair_spinaware_Kᵈ)(::Val{:buffer_uu})
	buffer_nw = Vector{ComplexF64}(undef, W.nw_up)
	buffer_mul_nw = Vector{ComplexF64}(undef, W.nw_up)
	buffer_np = Vector{ComplexF64}(undef, W.np_up)
	buffer_mul_np = Vector{ComplexF64}(undef, W.np_up)
	Up = Vector{ComplexF64}(undef, W.np_up_ij)
	return buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up
end
function (W::_U_pair_spinaware_Kᵈ)(::Val{:buffer_dd})
	buffer_nw = Vector{ComplexF64}(undef, W.nw_dn)
	buffer_mul_nw = Vector{ComplexF64}(undef, W.nw_dn)
	buffer_np = Vector{ComplexF64}(undef, W.np_dn)
	buffer_mul_np = Vector{ComplexF64}(undef, W.np_dn)
	Up = Vector{ComplexF64}(undef, W.np_dn_ij)
	return buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up
end
function (W::_U_pair_spinaware_Kᵈ)(::Val{:buffer_ud})
	buffer_nw_up = Vector{ComplexF64}(undef, W.nw_up)
	buffer_nw_dn = Vector{ComplexF64}(undef, W.nw_dn)
	buffer_mul_nw = Vector{ComplexF64}(undef, W.nw_dn)
	buffer_np_up = Vector{ComplexF64}(undef, W.np_up)
	buffer_np_dn = Vector{ComplexF64}(undef, W.np_dn)
	buffer_mul_np = Vector{ComplexF64}(undef, W.np_dn)
	Up_up = Vector{ComplexF64}(undef, W.np_up_ij)
	Up_dn = Vector{ComplexF64}(undef, W.np_dn_ij)
	return buffer_nw_up, buffer_nw_dn, buffer_mul_nw, buffer_np_up, buffer_np_dn, buffer_mul_np, Up_up, Up_dn
end
function (W::_U_pair_spinaware_Kᵈ)(::Val{:buffer_du})
	buffer_nw_up = Vector{ComplexF64}(undef, W.nw_up)
	buffer_nw_dn = Vector{ComplexF64}(undef, W.nw_dn)
	buffer_mul_nw = Vector{ComplexF64}(undef, W.nw_up)
	buffer_np_up = Vector{ComplexF64}(undef, W.np_up)
	buffer_np_dn = Vector{ComplexF64}(undef, W.np_dn)
	buffer_mul_np = Vector{ComplexF64}(undef, W.np_up)
	Up_up = Vector{ComplexF64}(undef, W.np_up_ij)
	Up_dn = Vector{ComplexF64}(undef, W.np_dn_ij)
	return buffer_nw_up, buffer_nw_dn, buffer_mul_nw, buffer_np_up, buffer_np_dn, buffer_mul_np, Up_up, Up_dn
end
function (W::_U_pair_spinaware_Kᵈ)(ψv′, ψc′, k′, ψv, ψc, k)
	buffer = W(Val(:buffer))
	return W(buffer..., ψv′, ψc′, k′, ψv, ψc, k)
end
function (W::_U_pair_spinaware_Kᵈ)(::Val{sym}, ψv′, ψc′, k′, ψv, ψc, k) where {sym}
	buffer = W(Val(Symbol(:buffer_, sym)))
	return W(Val(sym), buffer..., ψv′, ψc′, k′, ψv, ψc, k)
end
function (W::_U_pair_spinaware_Kᵈ)(
	buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up,
	ψv′, ψc′, k′, ψv, ψc, k,
)

	pair = W.pair
	pair_ij = W.pair_ij
	nw = W.nw

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
	Bw = view(Up, 1:nw)
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

	Aw = view(Up, 1:nw)
	Kᵈ = Aw ⋅ buffer_mul_nw
	Kᵈ += buffer_np ⋅ buffer_mul_np

	return -Kᵈ
end
function (W::_U_pair_spinaware_Kᵈ)(::Val{:uu},
	buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up,
	ψv′, ψc′, k′, ψv, ψc, k,
)

	pair = W.pair_up
	pair_ij = W.pair_up_ij
	nw = W.nw_up

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
	Bw = view(Up, 1:nw)
	Uk′k = view(W.Uk_uu,:,:,k′k)
	mul!(buffer_mul_nw, Uk′k, Bw)
	Pk′k = view(W.pairk_uu,:,:,k′k)
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

	Aw = view(Up, 1:nw)
	Kᵈ = Aw ⋅ buffer_mul_nw
	Kᵈ += buffer_np ⋅ buffer_mul_np

	return -Kᵈ
end
function (W::_U_pair_spinaware_Kᵈ)(::Val{:dd},
	buffer_nw, buffer_mul_nw, buffer_np, buffer_mul_np, Up,
	ψv′, ψc′, k′, ψv, ψc, k,
)

	pair = W.pair_dn
	pair_ij = W.pair_dn_ij
	nw = W.nw_dn

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
	Bw = view(Up, 1:nw)
	Uk′k = view(W.Uk_dd,:,:,k′k)
	mul!(buffer_mul_nw, Uk′k, Bw)
	Pk′k = view(W.pairk_dd,:,:,k′k)
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

	Aw = view(Up, 1:nw)
	Kᵈ = Aw ⋅ buffer_mul_nw
	Kᵈ += buffer_np ⋅ buffer_mul_np

	return -Kᵈ
end
function (W::_U_pair_spinaware_Kᵈ)(::Val{:ud},
	buffer_nw_up, buffer_nw_dn, buffer_mul_nw, buffer_np_up, buffer_np_dn, buffer_mul_np, Up_up, Up_dn,
	ψv′, ψc′, k′, ψv, ψc, k,
)

	pair = W.pair_up
	pair_ij = W.pair_up_ij
	nw = W.nw_up

	buffer_nw_up .= conj.(ψv)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up_up[idx] = ψv′[i] * buffer_nw_up[j]
	end
	ek_pR = view(W.ek_pR, :, k)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np_up[idx] = Up_up[ij] * ek_pR[R] # Bp
	end

	k′k = W.kgrid_minusmap[k′, k]
	Bw = view(Up_up, 1:nw)
	Uk′k = view(W.Uk_du,:,:,k′k)
	mul!(buffer_mul_nw, Uk′k, Bw)
	Pk′k = view(W.pairk_du,:,:,k′k)
	mul!(buffer_mul_np, Pk′k, buffer_np_up)

	pair = W.pair_dn
	pair_ij = W.pair_dn_ij
	nw = W.nw_dn

	buffer_nw_dn .= conj.(ψc)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up_dn[idx] = ψc′[i] * buffer_nw_dn[j]
	end
	ekq_pR = view(W.ekq_pR, :, k)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np_dn[idx] = Up_dn[ij] * ekq_pR[R] # Ap
	end

	Aw = view(Up_dn, 1:nw)
	Kᵈ = Aw ⋅ buffer_mul_nw
	Kᵈ += buffer_np_dn ⋅ buffer_mul_np

	return -Kᵈ
end
function (W::_U_pair_spinaware_Kᵈ)(::Val{:du},
	buffer_nw_up, buffer_nw_dn, buffer_mul_nw, buffer_np_up, buffer_np_dn, buffer_mul_np, Up_up, Up_dn,
	ψv′, ψc′, k′, ψv, ψc, k,
)

	pair = W.pair_dn
	pair_ij = W.pair_dn_ij
	nw = W.nw_dn

	buffer_nw_dn .= conj.(ψv)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up_dn[idx] = ψv′[i] * buffer_nw_dn[j]
	end
	ek_pR = view(W.ek_pR, :, k)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np_dn[idx] = Up_dn[ij] * ek_pR[R] # Bp
	end

	k′k = W.kgrid_minusmap[k′, k]
	Bw = view(Up_dn, 1:nw)
	Uk′k = view(W.Uk_ud,:,:,k′k)
	mul!(buffer_mul_nw, Uk′k, Bw)
	Pk′k = view(W.pairk_ud,:,:,k′k)
	mul!(buffer_mul_np, Pk′k, buffer_np_dn)

	pair = W.pair_up
	pair_ij = W.pair_up_ij
	nw = W.nw_up

	buffer_nw_up .= conj.(ψc)
	@inbounds @simd for idx in eachindex(pair_ij)
		(i, j) = pair_ij[idx]
		Up_up[idx] = ψc′[i] * buffer_nw_up[j]
	end
	ekq_pR = view(W.ekq_pR, :, k)
	@inbounds @simd for idx in eachindex(pair)
		(ij, R) = pair[idx]
		buffer_np_up[idx] = Up_up[ij] * ekq_pR[R] # Ap
	end

	Aw = view(Up_up, 1:nw)
	Kᵈ = Aw ⋅ buffer_mul_nw
	Kᵈ += buffer_np_up ⋅ buffer_mul_np

	return -Kᵈ
end
function _U_pair_spinaware_Kᵈ(kgrid::RedKgrid, U::DirectInterAction, pair_int::WanIntijklR; upindex::Vector{Int}, dnindex::Vector{Int})

	nk = length(kgrid)
	nw = numorb(U)
	nw == pair_int.nw || error("Mismatched U and pair interaction!")

	nw_up = length(upindex)
	nw_dn = length(dnindex)

	kgrid_Γ = RedKgrid(MonkhorstPack(kgrid.kgrid_size))
	kgrid_minusmap = kgridmap(kgrid_Γ, kgrid, kgrid, -)
	iΓ = kgrid_minusmap[1, 1]

	# pair_up, pair_dn
	pair = pair_int.pair
	pair_upindex = findall(pair) do (iw1, iw2, _)
		iw1 ∈ upindex && iw2 ∈ upindex
	end
	pair_dnindex = findall(pair) do (iw1, iw2, _)
		iw1 ∈ dnindex && iw2 ∈ dnindex
	end

	pair_rh = ReciprocalHoppings(pair_int)

	np = pair_int.np
	np_up = length(pair_upindex)
	np_dn = length(pair_dnindex)

	Uk = Array{ComplexF64}(undef, nw, nw, nk)
	pairk = Array{ComplexF64}(undef, np, np, nk)
	Threads.@threads for k in Base.OneTo(nk)
		Uk_view = view(Uk,:,:,k)
		U(Uk_view, kgrid_Γ[k], kgrid.kgrid_size; isΓ = (k == iΓ))
		Uk_view ./= nk
		pairk_view = view(pairk,:,:,k)
		pair_rh(pairk_view, kgrid_Γ[k])
		pairk_view ./= nk
	end

	Uk_uu = Uk[upindex, upindex, :]
	Uk_dd = Uk[dnindex, dnindex, :]
	Uk_ud = Uk[upindex, dnindex, :]
	Uk_du = Uk[dnindex, upindex, :]

	pairk_uu = pairk[pair_upindex, pair_upindex, :]
	pairk_dd = pairk[pair_dnindex, pair_dnindex, :]
	pairk_ud = pairk[pair_upindex, pair_dnindex, :]
	pairk_du = pairk[pair_dnindex, pair_upindex, :]

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

	# pair
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
	eq_pR = Vector{ComplexF64}(undef, pair_int.npR)
	ekq_pR = similar(ek_pR)

	return _U_pair_spinaware_Kᵈ(nw, nw_up, nw_dn, np, np_up, np_dn, np_ij, np_up_ij, np_dn_ij, nk,
		kgrid, kgrid_Γ, kgrid_minusmap, upindex, dnindex, pair_upindex, pair_dnindex, U, pair_rh,
		Uk, Uk_uu, Uk_dd, Uk_ud, Uk_du,
		pairk, pairk_uu, pairk_dd, pairk_ud, pairk_du,
		pair, pair_ij, pair_up, pair_up_ij, pair_dn, pair_dn_ij,
		pair_int.pR, ek_pR, eq_pR, ekq_pR)
end
