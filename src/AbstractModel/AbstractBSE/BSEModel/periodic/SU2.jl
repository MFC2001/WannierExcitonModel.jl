struct BaseBSESU2{
	TBT <: AbstractTightBindModel,
	KT <: AbstractKernalInterAction,
} <: AbstractBSE
	TB::TBT
	scissor::Float64
	kgrid::RedKgrid
	unitcell::Vector{ReducedCoordinates{Int}}
	vckmap::vckMap
	ijRmap::ijRMap
	bandk::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	bandkq::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	Kernal::KT
end
mutable struct BSESU2{T <: BaseBSESU2} <: AbstractBSE
	BSE::T
	isqgrid::Bool
end
function Base.getproperty(bse::BSESU2, sym::Symbol)
	if sym === :BSE
		return getfield(bse, :BSE)
	elseif sym === :isqgrid
		return getfield(bse, :isqgrid)
	else
		return getfield(bse.BSE, sym)
	end
end
function Base.setproperty!(bse::BSESU2, sym::Symbol, v)
	if sym === :isqgrid
		setfield!(bse, :isqgrid, v)
		if v
			@info "run Kernal(:initialize_qgrid)"
			bse.Kernal(Val(:initialize_qgrid))
		end
	else
		setfield!(bse, sym, v)
	end
end
function Base.show(io::IO, bse::BSESU2)
	print(io, "$(count(bse.TB.period)) dimensinal BSE model with $(numatom(bse.TB)) atoms and $(numorb(bse.TB)) orbitals.")
end
"""
"""
function BSESU2(TB::AbstractTightBindModel, Kernal::AbstractKernalInterAction;
	kgrid::RedKgrid, v, c, scissor::Real = 0, isqgrid::Bool = false)

	unitcell = gridindex(kgrid.kgrid_size)

	vckmap = vckMap(v, c, length(kgrid))
	ijRmap = ijRMap(numorb(TB), length(unitcell))

	bandk = BAND(kgrid, TB; vector = true)
	bandkq = similar(bandk)

	Kernal(Val(:initialize))

	bse = BaseBSESU2(TB, Float64(scissor), kgrid, unitcell, vckmap, ijRmap, bandk, bandkq, Kernal)
	if isqgrid
		bse.Kernal(Val(:initialize_qgrid))
	end
	return BSESU2(bse, isqgrid)
end
function (bse::BSESU2)(q::ReducedCoordinates)
	N = length(bse.vckmap)
	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)
	return bse(Htriplet, Hsinglet, q)
end
function (bse::BSESU2)(Htriplet, Hsinglet, q::ReducedCoordinates)
	_BSE_preprocess_q!(bse, q, Val(bse.isqgrid))
	return _BSE_Hamiltonian!(bse, Htriplet, Hsinglet)
end
function _BSE_preprocess_q!(bse::BSESU2, q, ::Val{true})
	if norm(q) < 1e-8
		q = _BSE_preprocess_eleband_q!(bse, q, Val(false))
		bse.Kernal(q)
	else
		kq_kindex, kΓq_kΓindex = _BSE_preprocess_eleband_qnotΓ!(bse, q, Val(true))
		bse.Kernal(kq_kindex, kΓq_kΓindex, q)
	end
	return nothing
end
function _BSE_preprocess_q!(bse::BSESU2, q, ::Val{false})
	q = _BSE_preprocess_eleband_q!(bse, q, Val(false))
	bse.Kernal(q)
	return nothing
end
function _BSE_preprocess_eleband_q!(bse::BSESU2, q, ::Val{true})
	if norm(q) < 1e-8
		q = _BSE_preprocess_eleband_q!(bse, q, Val(false))
	else
		_BSE_preprocess_eleband_qnotΓ!(bse, q, Val(true))
	end
	return q
end
function _BSE_preprocess_eleband_q!(bse::BSESU2, q, ::Val{false})
	# We can't calculate Γ directly.
	if norm(q) < 1e-8
		q = _BSE_shiftΓ(bse.TB.period)
	end
	bandkq = BAND(map(k -> k + q, bse.kgrid), bse.TB; vector = true)

	bse.bandkq .= bandkq
	return q
end
function _BSE_preprocess_eleband_qnotΓ!(bse::BSESU2, q, ::Val{true})
	kgrid = bse.kgrid
	kgrid_Γ = bse.Kernal.kgrid_Γ

	nk = length(kgrid)
	kq_kindex = Vector{Int}(undef, nk)
	kΓq_kΓindex = Vector{Int}(undef, nk)
	Threads.@threads for ki in Base.OneTo(nk)
		kq = kgrid[ki] + q
		kq_kindex[ki] = findfirst(k -> all(isinteger, k - kq), kgrid)
		kΓq = kgrid_Γ[ki] + q
		kΓq_kΓindex[ki] = findfirst(k -> all(isinteger, k - kΓq), kgrid_Γ)
	end

	bse.bandkq .= bse.bandk[kq_kindex]
	return kq_kindex, kΓq_kΓindex
end

function _BSE_Hamiltonian!(bse::BSESU2, Htriplet, Hsinglet)
	vckmap = bse.vckmap
	kernal = bse.Kernal
	bandk = bse.bandk
	bandkq = bse.bandkq
	scissor = bse.scissor

	N = length(vckmap)
	tasks = Vector{Task}(undef, Int(N * (N + 1) / 2))
	n = 0
	for j in 2:N, i in 1:j-1
		n += 1
		tasks[n] = Threads.@spawn begin
			(v′, c′, k′) = vckmap[i]
			(v, c, k) = vckmap[j]

			ψc′ = bandkq[k′].vectors[:, c′]
			ψc = bandkq[k].vectors[:, c]
			ψv′ = bandk[k′].vectors[:, v′]
			ψv = bandk[k].vectors[:, v]

			Kᵈ, Kˣ = kernal(k′, k, ψc′, ψv′, ψv, ψc)

			Htriplet[i, j] = Kᵈ
			Hsinglet[i, j] = Kᵈ + 2 * Kˣ
		end
	end
	for i in 1:N
		n += 1
		tasks[n] = Threads.@spawn begin
			(v, c, k) = vckmap[i]

			ψc = bandkq[k].vectors[:, c]
			ψv = bandk[k].vectors[:, v]

			Kᵈ, Kˣ = kernal(k, k, ψc, ψv, copy(ψv), copy(ψc))

			Δϵ = bandkq[k].values[c] - bandk[k].values[v] + scissor
			Htriplet[i, i] = Δϵ + Kᵈ
			Hsinglet[i, i] = Δϵ + Kᵈ + 2 * Kˣ
		end
	end
	wait.(tasks)

	return Hermitian(Htriplet, :U), Hermitian(Hsinglet, :U)
end

"""
	_uijR_ψvck(bse::BSESU2, η)

return an instance, which is used to extract periodic part of exciton bloch wave function.
"""
function _uijR_ψvck(bse::BSESU2, ηt, ηs)

	norb = length(bse.TB.orb_location)

	vckmap = bse.vckmap
	ijRmap = bse.ijRmap

	Nvck = length(vckmap)
	NijR = length(ijRmap)

	ekR = [cis(2π * (k ⋅ R)) for R in bse.unitcell, k in bse.kgrid]
	ekqR = similar(ekR)

	UU = Array{ComplexF64}(undef, norb, norb, Nvck)

	if ηt == ηs
		phase = _uijR_ψvck_phase(ηt, bse.unitcell, bse.TB.orb_location)
		BM = Matrix{ComplexF64}(undef, NijR, Nvck)
		return _uijR_ψvck_oneη(Nvck, NijR, vckmap, ijRmap, phase, ekR, ekqR, UU, BM)
	else
		phase_t = _uijR_ψvck_phase(ηt, bse.unitcell, bse.TB.orb_location)
		phase_s = _uijR_ψvck_phase(ηs, bse.unitcell, bse.TB.orb_location)
		BMt = Matrix{ComplexF64}(undef, NijR, Nvck)
		BMs = Matrix{ComplexF64}(undef, NijR, Nvck)
		return _uijR_ψvck_twoη(Nvck, NijR, vckmap, ijRmap, phase_t, phase_s, ekR, ekqR, UU, BMt, BMs)
	end
end
