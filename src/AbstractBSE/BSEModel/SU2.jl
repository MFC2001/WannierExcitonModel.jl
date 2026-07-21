struct BSESU2{
	TBT <: AbstractTightBindModel,
	KT <: KernelInterAction,
} <: AbstractBSE
	TB::TBT
	scissor::Float64
	kgrid::RedKgrid
	unitcell::Vector{ReducedCoordinates{Int}}
	vckmap::vckMap
	ijRmap::ijRMap
	bandk::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	bandkq::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	kernel::KT
	Γdata::_BSE_Γdata
	H_upidx::Vector{NTuple{2, Int}}
end
function Base.show(io::IO, bse::BSESU2)
	print(io, "$(count(bse.TB.period)) dimensinal BSE model with $(numatom(bse.TB)) atoms and $(numorb(bse.TB)) orbitals.")
end
"""
"""
function BSE(::Val{:SU2}, TB::AbstractTightBindModel, kernel::KernelInterAction; v, c, scissor::Real = 0)

	kgrid = kernel.kgrid
	unitcell = gridindex(Val(:WignerSeitz), TB.lattice, kgrid.kgrid_size)

	@assert length(v) > 0 "Empty `v` is unaccepted!"
	@assert length(c) > 0 "Empty `c` is unaccepted!"
	vckmap = vckMap(v, c, length(kgrid))
	ijRmap = ijRMap(numorb(TB), length(unitcell))

	bandk = BAND(kgrid, TB; vector = true)
	phase_normalize!(bandk, :real_sum)
	bandkq = similar(bandk)

	Γdata = _BSE_Γdata(TB, kgrid, bandk, vckmap)

	N = length(vckmap)
	H_upidx = Vector{NTuple{2, Int}}(undef, Int(N * (N + 1) / 2))
	n = 0
	for j in 1:N, i in 1:j
		n += 1
		H_upidx[n] = (i, j)
	end

	return BSESU2(TB, Float64(scissor), kgrid, unitcell, vckmap, ijRmap, bandk, bandkq, kernel, Γdata, H_upidx)
end
function (bse::BSESU2)(::Val{:dimension})
	return length(bse.vckmap)
end
# isΓ表示是否精确处理Γ点，这时输入的q表示方向，注意是分数坐标。
function (bse::BSESU2)(q::ReducedCoordinates; isqgrid::Bool = false, isΓ::Bool = false)
	N = length(bse.vckmap)
	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)
	return bse(Htriplet, Hsinglet, q; isqgrid, isΓ)
end
function (bse::BSESU2)(Htriplet, Hsinglet, q::ReducedCoordinates; isqgrid::Bool = false, isΓ::Bool = false)
	(q, isΓ) = _BSE_preprocess_Γ(q, isΓ, bse.TB.period)
	_BSE_preprocess_q!(bse, q; isqgrid, isΓ)
	return _BSE_Hamiltonian!(Htriplet, Hsinglet, bse, q; isΓ)
end
"""
if isΓ, don't use q.
核心是kernel，其设计中不计入方向依赖项。
"""
function _BSE_preprocess_q!(bse::BSESU2, q::ReducedCoordinates; isqgrid::Bool, isΓ::Bool)
	if isΓ
		bse.bandkq .= bse.bandk
		bse.kernel(ReducedCoordinates(0, 0, 0); isΓ = true)
	else
		if isqgrid
			kgrid = bse.kgrid
			bandk = bse.bandk
			bandkq = bse.bandkq
			Threads.@threads for ik in eachindex(kgrid)
				kq = kgrid[ik] + q
				kq_kindex = findfirst(k -> all(isinteger, k - kq), kgrid)
				bandkq[ik] = bandk[kq_kindex]
			end
		else
			bandkq = BAND(map(k -> k + q, bse.kgrid), bse.TB; vector = true)
			phase_normalize!(bandkq, :real_sum)
			bse.bandkq .= bandkq
		end
		bse.kernel(q; isΓ = false)
	end
	return nothing
end
function _BSE_Hamiltonian!(Htriplet, Hsinglet, bse::BSESU2, q::ReducedCoordinates; isΓ)
	vckmap = bse.vckmap
	kernel = bse.kernel
	bandk = bse.bandk
	bandkq = bse.bandkq
	scissor = bse.scissor

	LinearAlgebra.BLAS.set_num_threads(1)

	# buffers
	nt = Threads.nthreads()
	nchunk = nt * 2
	buffer_chunk = [kernel(Val(:buffer)) for _ in Base.OneTo(nchunk)]
	H_idx = _split_tasks(bse.H_upidx, nchunk)

	Threads.@threads for ichunk in Base.OneTo(nchunk)
		(buffer_Kᵈ, buffer_Kˣ) = buffer_chunk[ichunk]
		@inbounds for (i, j) in H_idx[ichunk]
			(v′, c′, k′) = vckmap[i]
			(v, c, k) = vckmap[j]
			ψv′ = view(bandk[k′].vectors, :, v′)
			ψc′ = view(bandkq[k′].vectors, :, c′)
			ψv = view(bandk[k].vectors, :, v)
			ψc = view(bandkq[k].vectors, :, c)
			Kᵈ, Kˣ = kernel(buffer_Kᵈ, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)

			if isfinite(Kᵈ) && isfinite(Kˣ)
				if i == j
					Δϵ = bandkq[k].values[c] - bandk[k].values[v] + scissor
					Htriplet[i, i] = Δϵ + Kᵈ
					Hsinglet[i, i] = Δϵ + Kᵈ + 2 * Kˣ
				else
					Htriplet[i, j] = Kᵈ
					Hsinglet[i, j] = Kᵈ + 2 * Kˣ
				end
			else
				throw(DomainError(
					(i = i, j = j, k′ = k′, k = k, Kᵈ = Kᵈ, Kˣ = Kˣ),
					"BSE matrix element at isΓ=$isΓ, q=$q, (i,j)=($i,$j) has NaN/Inf: Kᵈ=$Kᵈ, Kˣ=$Kˣ",
				))
			end
		end
	end

	if isΓ && count(bse.TB.period) == 3
		_BSE_KˣΓ(bse, Htriplet, Hsinglet, H_idx, q)
	end

	LinearAlgebra.BLAS.set_num_threads(nt)

	return Hermitian(Htriplet, :U), Hermitian(Hsinglet, :U)
end
function _BSE_KˣΓ(bse::BSESU2, Htriplet, Hsinglet, H_idx, q)
	qcar = bse.Γdata.rlattice * q
	qnorm = norm(qcar)
	q̂x = qcar[1] / qnorm
	q̂y = qcar[2] / qnorm
	q̂z = qcar[3] / qnorm
	upu = bse.Γdata.upu
	Threads.@threads for idxs in H_idx
		for (i, j) in idxs
			Hsinglet[i, j] += 2 * conj(upu[1, i] * q̂x + upu[2, i] * q̂y + upu[3, i] * q̂z) *
							  (upu[1, j] * q̂x + upu[2, j] * q̂y + upu[3, j] * q̂z) * bse.Γdata.CoulombScaledivnk
		end
	end
	return Htriplet, Hsinglet
end
