function _mmn_BSE_uband(qgrid, bse, uband, η, bandindex, nnkpts)

	kgrid = bse.kgrid
	vckmap = bse.vckmap
	(nv, nc, nk) = size(vckmap)

	nq = length(BSEband)
	nband = length(bandindex)
	norb_ele = length(bse.TB.orb_location)
	nqb = size(nnkpts, 2)


	k_plus_q_kindex = Matrix{Int}(undef, nk, nqb)
	k_plus_q_plus_b2_kindex = Matrix{Int}(undef, nk, nqb)
	k_minus_b2_kindex = Matrix{Int}(undef, nk, nqb)
	b2 = Vector{Vec3{Rational{Int}}}(undef, nqb)
	Threads.@threads for i in 1:nqb
		q_qindex = nnkpts[1, i]
		q = qgrid.kdirect[q_qindex]
		q_plus_b_qindex = nnkpts[2, i]
		q_plus_b = qgrid.kdirect[q_plus_b_qindex] + nnkpts[3:5, i]
		b2[i] = (q_plus_b - q) .// 2
		b2_p = [0.0, 0.0, 0.0]
		b2_p[Tp] .= b2[Tp]
		for (ki, k) in enumerate(kgrid)
			# k = kgrid[ki]

			k_plus_q = k + q
			k_plus_q_kindex[ki, i] = findfirst(k -> all(isinteger, k - k_plus_q), kgrid)

			k_plus_q_plus_b2 = k + q_plus_b - b2_p
			k_plus_q_plus_b2_kindex[ki, i] = findfirst(k -> all(isinteger, k - k_plus_q_plus_b2), kgrid)

			k_minus_b2 = k - b2_p
			k_minus_b2_kindex[ki, i] = findfirst(k -> all(isinteger, k - k_minus_b2), kgrid)
		end
	end

	#Calculate Mc and Mv
	eb2orb = Matrix{ComplexF64}(undef, norb_ele, nqb)
	Threads.@threads for i in 1:nqb
		eb2orb[:, i] = map(x -> cis(-2π * (b2[i] ⋅ x)), bse.TB.orb_location)
	end
	band = bse.bandk
	#Mc
	tasks = Array{Task}(undef, nc, nc, nk, nqb)
	for i in 1:nqb, k in 1:nk, c′ in 1:nc, c in 1:nc
		tasks[c, c′, k, i] = Threads.@spawn sum(Base.OneTo(norb_ele)) do iorb
			eb2orb[iorb, i] * conj(band[k_plus_q_kindex[k, i]].vectors[iorb, vckmap.idx2c[c]]) * band[k_plus_q_plus_b2_kindex[k, i]].vectors[iorb, vckmap.idx2c[c′]]
		end
	end
	Mc = fetch.(tasks)
	#Mv
	tasks = Array{Task}(undef, nv, nv, nk, nqb)
	for i in 1:nqb, k in 1:nk, v′ in 1:nv, v in 1:nv
		tasks[v, v′, k, i] = Threads.@spawn sum(Base.OneTo(norb_ele)) do iorb
			eb2orb[iorb, i] * conj(band[k_minus_b2_kindex[k, i]].vectors[iorb, vckmap.idx2v[v′]]) * band[k].vectors[iorb, vckmap.idx2v[v]]
		end
	end
	Mv = fetch.(tasks)

	#Calculate M
	nq = length(BSEband)
	Avck = Array{ComplexF64}(undef, nv, nc, nk, nband, nq)
	for q in 1:nq, α in 1:nband, k in 1:nk, c in 1:nc, v in 1:nv
		Avck[v, c, k, α, q] = BSEband[q].vectors[vckmap.vck2idx[v, c, k], bandindex[α]]
	end

	vvcck = Tuple.(CartesianIndices((nv, nv, nc, nc, nk)))




	norb = length(bse.TB.orb_location)
	nR = length(bse.unitcell)

	δ = 1 - η
	orb₀ = map(CartesianIndices((norb, norb))) do I
		(i, j) = Tuple(I)
		(δ * bse.TB.orb_location[i] + η * bse.TB.orb_location[j])
	end

	eqorb₀ = Matrix{ComplexF64}(undef, norb, norb)
	eqR = Vector{ComplexF64}(undef, nR)

	tasks = Array{Task}(undef, nband, nband, nqb)
	for i in 1:nqb, β in 1:nband, α in 1:nband
		tasks[α, β, i] = Threads.@spawn begin
			G = ReducedCoordinates(nnkpts[3:5, i])
			if iszero(G)
				M = uband[q].vectors[:, α] ⋅ uband[qb].vectors[:, β]
			else
				M = uband[q].vectors[:, α] ⋅ uband[qb].vectors[:, β]
			end
			M
		end
	end
	M = fetch.(tasks)

	return M
end
