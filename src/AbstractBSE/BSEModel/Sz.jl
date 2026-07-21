struct BSESz{
	TBT_up <: AbstractTightBindModel,
	TBT_dn <: AbstractTightBindModel,
	KT <: KernelInterAction,
} <: AbstractBSE
	TB_up::TBT_up
	TB_dn::TBT_dn
	scissor_c_up::Float64
	scissor_v_up::Float64
	scissor_c_dn::Float64
	scissor_v_dn::Float64
	kgrid::RedKgrid
	unitcell::Vector{ReducedCoordinates{Int}}
	vckmap_uu::vckMap
	vckmap_dd::vckMap
	vckmap_ud::vckMap
	vckmap_du::vckMap
	ijRmap_uu::ijRMap
	ijRmap_dd::ijRMap
	ijRmap_ud::ijRMap
	ijRmap_du::ijRMap
	bandk_up::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	bandk_dn::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	bandkq_up::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	bandkq_dn::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	kernel::KT
	Γdata_up::_BSE_Γdata
	Γdata_dn::_BSE_Γdata
	Huu_upidx::Vector{NTuple{2, Int}}
	Hdd_upidx::Vector{NTuple{2, Int}}
	Huudd_idx::Vector{NTuple{2, Int}}
	Hud_upidx::Vector{NTuple{2, Int}}
	Hdu_upidx::Vector{NTuple{2, Int}}
end
function Base.show(io::IO, bse::BSESz)
	print(io, "$(count(bse.TB_up.period)) dimensinal BSE model with $(numorb(bse.TB_up)) up orbitals and $(numorb(bse.TB_dn)) down orbitals.")
end
"""
"""
function BSE(::Val{:Sz}, TB_up::AbstractTightBindModel, TB_dn::AbstractTightBindModel, kernel::KernelInterAction;
	v_up = Int[], c_up = Int[], v_dn = Int[], c_dn = Int[],
	scissor::Real = 0, scissor_up::Real = scissor, scissor_dn::Real = scissor,
	scissor_c_up::Real = scissor_up, scissor_c_dn::Real = scissor_dn, scissor_v_up::Real = 0, scissor_v_dn::Real = 0)

	@info "If the interaction parameters provided are defined on the basis of wannier functions of **different spin**, please set Kernel(...; interaction = :spinaware)."

	kgrid = kernel.kgrid
	unitcell = gridindex(Val(:WignerSeitz), TB_up.lattice, kgrid.kgrid_size)

	nk = length(kgrid)
	vckmap_uu = vckMap(v_up, c_up, nk)
	vckmap_dd = vckMap(v_dn, c_dn, nk)
	vckmap_ud = vckMap(v_up, c_dn, nk)
	vckmap_du = vckMap(v_dn, c_up, nk)

	nw_up = numorb(TB_up)
	nw_dn = numorb(TB_dn)
	nR = length(unitcell)
	ijRmap_uu = ijRMap(nw_up, nw_up, nR)
	ijRmap_dd = ijRMap(nw_dn, nw_dn, nR)
	ijRmap_ud = ijRMap(nw_up, nw_dn, nR)
	ijRmap_du = ijRMap(nw_dn, nw_up, nR)

	bandk_up = BAND(kgrid, TB_up; vector = true)
	bandk_dn = BAND(kgrid, TB_dn; vector = true)
	phase_normalize!(bandk_up, :real_sum)
	phase_normalize!(bandk_dn, :real_sum)
	bandkq_up = similar(bandk_up)
	bandkq_dn = similar(bandk_dn)

	Γdata_up = _BSE_Γdata(TB_up, kgrid, bandk_up, vckmap_uu)
	Γdata_dn = _BSE_Γdata(TB_dn, kgrid, bandk_dn, vckmap_dd)

	Nuu = length(vckmap_uu)
	Ndd = length(vckmap_dd)
	Nud = length(vckmap_ud)
	Ndu = length(vckmap_du)

	Huu_upidx = Vector{NTuple{2, Int}}(undef, Int(Nuu * (Nuu + 1) / 2))
	Hdd_upidx = Vector{NTuple{2, Int}}(undef, Int(Ndd * (Ndd + 1) / 2))
	Huudd_idx = Vector{NTuple{2, Int}}(undef, Nuu * Ndd)
	Hud_upidx = Vector{NTuple{2, Int}}(undef, Int(Nud * (Nud + 1) / 2))
	Hdu_upidx = Vector{NTuple{2, Int}}(undef, Int(Ndu * (Ndu + 1) / 2))

	n = 0
	for j in 1:Nuu, i in 1:j
		n += 1
		Huu_upidx[n] = (i, j)
	end
	n = 0
	for j in (Nuu+1):(Nuu+Ndd), i in (Nuu+1):j
		n += 1
		Hdd_upidx[n] = (i, j)
	end
	n = 0
	for j in (Nuu+1):(Nuu+Ndd), i in 1:Nuu
		n += 1
		Huudd_idx[n] = (i, j)
	end
	n = 0
	for j in 1:Nud, i in 1:j
		n += 1
		Hud_upidx[n] = (i, j)
	end
	n = 0
	for j in 1:Ndu, i in 1:j
		n += 1
		Hdu_upidx[n] = (i, j)
	end

	return BSESz(TB_up, TB_dn,
		Float64(scissor_c_up), Float64(scissor_v_up),
		Float64(scissor_c_dn), Float64(scissor_v_dn),
		kgrid, unitcell,
		vckmap_uu, vckmap_dd, vckmap_ud, vckmap_du,
		ijRmap_uu, ijRmap_dd, ijRmap_ud, ijRmap_du,
		bandk_up, bandk_dn, bandkq_up, bandkq_dn,
		kernel, Γdata_up, Γdata_dn,
		Huu_upidx, Hdd_upidx, Huudd_idx, Hud_upidx, Hdu_upidx,
	)
end
function (bse::BSESz)(::Val{:dimension})
	return length(bse.vckmap_uu) + length(bse.vckmap_dd), length(bse.vckmap_du), length(bse.vckmap_ud)
end
function (bse::BSESz)(q::ReducedCoordinates; isqgrid::Bool = false, isΓ::Bool = false)
	N0, Np1, Nn1 = bse(Val(:dimension))
	H0 = Matrix{ComplexF64}(undef, N0, N0)
	Hp1 = Matrix{ComplexF64}(undef, Np1, Np1)
	Hn1 = Matrix{ComplexF64}(undef, Nn1, Nn1)
	return bse(H0, Hp1, Hn1, q; isqgrid, isΓ)
end
function (bse::BSESz)(H0, Hp1, Hn1, q::ReducedCoordinates; isqgrid::Bool = false, isΓ::Bool = false)
	(q, isΓ) = _BSE_preprocess_Γ(q, isΓ, bse.TB_up.period)
	_BSE_preprocess_q!(bse, q; isqgrid, isΓ)
	return _BSE_Hamiltonian!(H0, Hp1, Hn1, bse, q; isΓ)
end
function _BSE_preprocess_q!(bse::BSESz, q::ReducedCoordinates; isqgrid::Bool, isΓ::Bool)
	if isΓ
		bse.bandkq_up .= bse.bandk_up
		bse.bandkq_dn .= bse.bandk_dn
		bse.kernel(ReducedCoordinates(0, 0, 0); isΓ = true)
	else
		if isqgrid
			kgrid = bse.kgrid
			bandk_up = bse.bandk_up
			bandk_dn = bse.bandk_dn
			bandkq_up = bse.bandkq_up
			bandkq_dn = bse.bandkq_dn
			Threads.@threads for ik in eachindex(kgrid)
				kq = kgrid[ik] + q
				kq_kindex = findfirst(k -> all(isinteger, k - kq), kgrid)
				bandkq_up[ik] = bandk_up[kq_kindex]
				bandkq_dn[ik] = bandk_dn[kq_kindex]
			end
		else
			kq = map(k -> k + q, bse.kgrid)
			bandkq_up = BAND(kq, bse.TB_up; vector = true)
			bandkq_dn = BAND(kq, bse.TB_dn; vector = true)
			phase_normalize!(bandkq_up, :real_sum)
			phase_normalize!(bandkq_dn, :real_sum)
			bse.bandkq_up .= bandkq_up
			bse.bandkq_dn .= bandkq_dn
		end
		bse.kernel(q; isΓ = false)
	end
	return nothing
end
function _BSE_Hamiltonian!(H0, Hp1, Hn1, bse::BSESz, q::ReducedCoordinates; isΓ)
	kernel = bse.kernel
	bandk_up = bse.bandk_up
	bandk_dn = bse.bandk_dn
	bandkq_up = bse.bandkq_up
	bandkq_dn = bse.bandkq_dn
	scissor_c_up = bse.scissor_c_up
	scissor_v_up = bse.scissor_v_up
	scissor_c_dn = bse.scissor_c_dn
	scissor_v_dn = bse.scissor_v_dn

	LinearAlgebra.BLAS.set_num_threads(1)

	# buffers
	nt = Threads.nthreads()
	nchunk = nt

	vckmap_uu = bse.vckmap_uu
	vckmap_dd = bse.vckmap_dd
	N_uu = length(vckmap_uu)

	if N_uu > 0
		H_uu_idx = _split_tasks(bse.Huu_upidx, nchunk)
		buffer_chunk_uu = [kernel(Val(:buffer_uu)) for _ in Base.OneTo(nchunk)]
		scissor_uu = scissor_c_up - scissor_v_up
		Threads.@threads for ichunk in Base.OneTo(nchunk)
			(buffer_Kᵈ, buffer_Kˣ) = buffer_chunk_uu[ichunk]
			@inbounds for (i, j) in H_uu_idx[ichunk]
				(v′, c′, k′) = vckmap_uu[i]
				(v, c, k) = vckmap_uu[j]
				ψv′ = view(bandk_up[k′].vectors, :, v′)
				ψc′ = view(bandkq_up[k′].vectors, :, c′)
				ψv = view(bandk_up[k].vectors, :, v)
				ψc = view(bandkq_up[k].vectors, :, c)
				Kᵈ, Kˣ = kernel(Val(:uu), buffer_Kᵈ, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)

				if isfinite(Kᵈ) && isfinite(Kˣ)
					if i == j
						Δϵ = bandkq_up[k].values[c] - bandk_up[k].values[v] + scissor_uu
						H0[i, i] = Δϵ + Kᵈ + Kˣ
					else
						H0[i, j] = Kᵈ + Kˣ
					end
				else
					throw(DomainError(
						(i = i, j = j, k′ = k′, k = k, Kᵈ = Kᵈ, Kˣ = Kˣ),
						"BSE matrix element at isΓ=$isΓ, q=$q, (i,j)=($i,$j) has NaN/Inf: Kᵈ=$Kᵈ, Kˣ=$Kˣ",
					))
				end
			end
		end
		if isΓ && count(bse.TB_up.period) == 3
			_BSE_KˣΓ(bse, Val(:uu), H0, H_uu_idx, q)
		end
	end

	if length(vckmap_dd) > 0
		H_dd_idx = _split_tasks(bse.Hdd_upidx, nchunk)
		buffer_chunk_dd = [kernel(Val(:buffer_dd)) for _ in Base.OneTo(nchunk)]
		scissor_dd = scissor_c_dn - scissor_v_dn
		Threads.@threads for ichunk in Base.OneTo(nchunk)
			(buffer_Kᵈ, buffer_Kˣ) = buffer_chunk_dd[ichunk]
			@inbounds for (i, j) in H_dd_idx[ichunk]
				(v′, c′, k′) = vckmap_dd[i-N_uu]
				(v, c, k) = vckmap_dd[j-N_uu]
				ψv′ = view(bandk_dn[k′].vectors, :, v′)
				ψc′ = view(bandkq_dn[k′].vectors, :, c′)
				ψv = view(bandk_dn[k].vectors, :, v)
				ψc = view(bandkq_dn[k].vectors, :, c)
				Kᵈ, Kˣ = kernel(Val(:dd), buffer_Kᵈ, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)

				if isfinite(Kᵈ) && isfinite(Kˣ)
					if i == j
						Δϵ = bandkq_dn[k].values[c] - bandk_dn[k].values[v] + scissor_dd
						H0[i, i] = Δϵ + Kᵈ + Kˣ
					else
						H0[i, j] = Kᵈ + Kˣ
					end
				else
					throw(DomainError(
						(i = i, j = j, k′ = k′, k = k, Kᵈ = Kᵈ, Kˣ = Kˣ),
						"BSE matrix element at isΓ=$isΓ, q=$q, (i,j)=($i,$j) has NaN/Inf: Kᵈ=$Kᵈ, Kˣ=$Kˣ",
					))
				end
			end
		end
		if isΓ && count(bse.TB_up.period) == 3
			_BSE_KˣΓ(bse, Val(:dd), H0, H_dd_idx, q)
		end
	end

	if length(vckmap_uu) > 0 && length(vckmap_dd) > 0
		H_uudd_idx = _split_tasks(bse.Huudd_idx, nchunk)
		buffer_chunk_uudd = [kernel(Val(:buffer_uudd)) for _ in Base.OneTo(nchunk)]
		Threads.@threads for ichunk in Base.OneTo(nchunk)
			buffer_Kˣ = buffer_chunk_uudd[ichunk]
			@inbounds for (i, j) in H_uudd_idx[ichunk]
				(v′, c′, k′) = vckmap_uu[i]
				(v, c, k) = vckmap_dd[j-N_uu]
				ψv′ = view(bandk_up[k′].vectors, :, v′)
				ψc′ = view(bandkq_up[k′].vectors, :, c′)
				ψv = view(bandk_dn[k].vectors, :, v)
				ψc = view(bandkq_dn[k].vectors, :, c)
				Kˣ = kernel(Val(:uudd), buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)

				if isfinite(Kˣ)
					H0[i, j] = Kˣ
				else
					throw(DomainError(
						(i = i, j = j, k′ = k′, k = k, Kˣ = Kˣ),
						"BSE matrix element at isΓ=$isΓ, q=$q, (i,j)=($i,$j) has NaN/Inf: Kᵈ=$Kᵈ, Kˣ=$Kˣ",
					))
				end
			end
		end
		if isΓ && count(bse.TB_up.period) == 3
			_BSE_KˣΓ(bse, Val(:uudd), H0, H_uudd_idx, q)
		end
	end

	vckmap_ud = bse.vckmap_ud
	if length(vckmap_ud) > 0
		H_ud_idx = _split_tasks(bse.Hud_upidx, nchunk)
		buffer_chunk_ud = [kernel(Val(:buffer_ud)) for _ in Base.OneTo(nchunk)]
		scissor_ud = scissor_c_dn - scissor_v_up
		Threads.@threads for ichunk in Base.OneTo(nchunk)
			buffer_Kᵈ = buffer_chunk_ud[ichunk]
			@inbounds for (i, j) in H_ud_idx[ichunk]
				(v′, c′, k′) = vckmap_ud[i]
				(v, c, k) = vckmap_ud[j]
				ψv′ = view(bandk_up[k′].vectors, :, v′)
				ψc′ = view(bandkq_dn[k′].vectors, :, c′)
				ψv = view(bandk_up[k].vectors, :, v)
				ψc = view(bandkq_dn[k].vectors, :, c)
				Kᵈ = kernel(Val(:ud), buffer_Kᵈ, ψv′, ψc′, k′, ψv, ψc, k)

				if isfinite(Kᵈ)
					if i == j
						Δϵ = bandkq_dn[k].values[c] - bandk_up[k].values[v] + scissor_ud
						Hn1[i, i] = Δϵ + Kᵈ
					else
						Hn1[i, j] = Kᵈ
					end
				else
					throw(DomainError(
						(i = i, j = j, k′ = k′, k = k, Kᵈ = Kᵈ),
						"BSE matrix element at isΓ=$isΓ, q=$q, (i,j)=($i,$j) has NaN/Inf: Kᵈ=$Kᵈ, Kˣ=$Kˣ",
					))
				end
			end
		end
	end


	vckmap_du = bse.vckmap_du
	if length(vckmap_du) > 0
		H_du_idx = _split_tasks(bse.Hdu_upidx, nchunk)
		buffer_chunk_du = [kernel(Val(:buffer_du)) for _ in Base.OneTo(nchunk)]
		scissor_du = scissor_c_up - scissor_v_dn
		Threads.@threads for ichunk in Base.OneTo(nchunk)
			buffer_Kᵈ = buffer_chunk_du[ichunk]
			@inbounds for (i, j) in H_du_idx[ichunk]
				(v′, c′, k′) = vckmap_du[i]
				(v, c, k) = vckmap_du[j]
				ψv′ = view(bandk_dn[k′].vectors, :, v′)
				ψc′ = view(bandkq_up[k′].vectors, :, c′)
				ψv = view(bandk_dn[k].vectors, :, v)
				ψc = view(bandkq_up[k].vectors, :, c)
				Kᵈ = kernel(Val(:du), buffer_Kᵈ, ψv′, ψc′, k′, ψv, ψc, k)

				if isfinite(Kᵈ)
					if i == j
						Δϵ = bandkq_up[k].values[c] - bandk_dn[k].values[v] + scissor_du
						Hp1[i, i] = Δϵ + Kᵈ
					else
						Hp1[i, j] = Kᵈ
					end
				else
					throw(DomainError(
						(i = i, j = j, k′ = k′, k = k, Kᵈ = Kᵈ),
						"BSE matrix element at isΓ=$isΓ, q=$q, (i,j)=($i,$j) has NaN/Inf: Kᵈ=$Kᵈ, Kˣ=$Kˣ",
					))
				end
			end
		end
	end

	LinearAlgebra.BLAS.set_num_threads(nt)

	return Hermitian(H0, :U), Hermitian(Hp1, :U), Hermitian(Hn1, :U)
end
function _BSE_KˣΓ(bse::BSESz, ::Val{:uu}, H0, H_idx, q)
	qcar = bse.Γdata_up.rlattice * q
	qnorm = norm(qcar)
	q̂x = qcar[1] / qnorm
	q̂y = qcar[2] / qnorm
	q̂z = qcar[3] / qnorm
	upu_uu = bse.Γdata_up.upu
	Threads.@threads for idxs in H_idx
		for (i, j) in idxs
			H0[i, j] += conj(upu_uu[1, i] * q̂x + upu_uu[2, i] * q̂y + upu_uu[3, i] * q̂z) *
						(upu_uu[1, j] * q̂x + upu_uu[2, j] * q̂y + upu_uu[3, j] * q̂z) * bse.Γdata_up.CoulombScaledivnk
		end
	end
	return H0
end
function _BSE_KˣΓ(bse::BSESz, ::Val{:dd}, H0, H_idx, q)
	qcar = bse.Γdata_up.rlattice * q
	qnorm = norm(qcar)
	q̂x = qcar[1] / qnorm
	q̂y = qcar[2] / qnorm
	q̂z = qcar[3] / qnorm
	upu_dd = bse.Γdata_dn.upu
	Threads.@threads for idxs in H_idx
		for (i, j) in idxs
			H0[i, j] += conj(upu_dd[1, i] * q̂x + upu_dd[2, i] * q̂y + upu_dd[3, i] * q̂z) *
						(upu_dd[1, j] * q̂x + upu_dd[2, j] * q̂y + upu_dd[3, j] * q̂z) * bse.Γdata_up.CoulombScaledivnk
		end
	end
	return H0
end
function _BSE_KˣΓ(bse::BSESz, ::Val{:uudd}, H0, H_idx, q)
	qcar = bse.Γdata_up.rlattice * q
	qnorm = norm(qcar)
	q̂x = qcar[1] / qnorm
	q̂y = qcar[2] / qnorm
	q̂z = qcar[3] / qnorm
	upu_uu = bse.Γdata_up.upu
	upu_dd = bse.Γdata_dn.upu
	Threads.@threads for idxs in H_idx
		for (i, j) in idxs
			H0[i, j] += conj(upu_uu[1, i] * q̂x + upu_uu[2, i] * q̂y + upu_uu[3, i] * q̂z) *
						(upu_dd[1, j] * q̂x + upu_dd[2, j] * q̂y + upu_dd[3, j] * q̂z) * bse.Γdata_up.CoulombScaledivnk
		end
	end
	return H0
end
function _BSE_KˣΓ(bse::BSESz, ::Val{:dduu}, H0, H_idx, q)
	qcar = bse.Γdata_up.rlattice * q
	qnorm = norm(qcar)
	q̂x = qcar[1] / qnorm
	q̂y = qcar[2] / qnorm
	q̂z = qcar[3] / qnorm
	upu_uu = bse.Γdata_up.upu
	upu_dd = bse.Γdata_dn.upu
	Threads.@threads for idxs in H_idx
		for (i, j) in idxs
			H0[i, j] += conj(upu_dd[1, i] * q̂x + upu_dd[2, i] * q̂y + upu_dd[3, i] * q̂z) *
						(upu_uu[1, j] * q̂x + upu_uu[2, j] * q̂y + upu_uu[3, j] * q̂z) * bse.Γdata_up.CoulombScaledivnk
		end
	end
	return H0
end
