struct BaseBSEgeneral{
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
mutable struct BSEgeneral{T <: BaseBSEgeneral} <: AbstractBSE
	BSE::T
	isqgrid::Bool
end
function Base.getproperty(bse::BSEgeneral, sym::Symbol)
	if sym === :BSE
		return getfield(bse, :BSE)
	elseif sym === :isqgrid
		return getfield(bse, :isqgrid)
	else
		return getfield(bse.BSE, sym)
	end
end
function Base.setproperty!(bse::BSEgeneral, sym::Symbol, v)
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
function Base.show(io::IO, bse::BSEgeneral)
	print(io, "$(count(bse.TB.period)) dimensinal BSE model with $(numatom(bse.TB)) atoms and $(numorb(bse.TB)) orbitals.")
end
"""
"""
function BSEgeneral(TB::AbstractTightBindModel, Kernal::AbstractKernalInterAction;
	kgrid::RedKgrid, v, c, scissor::Real = 0, isqgrid::Bool = false)

	unitcell = gridindex(kgrid.kgrid_size)

	vckmap = vckMap(v, c, length(kgrid))
	ijRmap = ijRMap(numorb(TB), length(unitcell))

	bandk = BAND(kgrid, TB; vector = true)
	bandkq = similar(bandk)

	Kernal(Val(:initialize))

	bse = BaseBSEgeneral(TB, Float64(scissor), kgrid, unitcell, vckmap, ijRmap, bandk, bandkq, Kernal)
	if isqgrid
		bse.Kernal(Val(:initialize_qgrid))
	end
	return BSEgeneral(bse, isqgrid)
end
function (bse::BSEgeneral)(q::ReducedCoordinates)
	N = length(bse.vckmap)
	H = Matrix{ComplexF64}(undef, N, N)
	return bse(H, q)
end
function (bse::BSEgeneral)(H, q::ReducedCoordinates)
	_BSE_preprocess_q!(bse, q, Val(bse.isqgrid))
	return _BSE_Hamiltonian!(bse, H)
end
function _BSE_preprocess_q!(bse::BSEgeneral, q, ::Val{true})
	if norm(q) < 1e-8
		q = _BSE_preprocess_eleband_q!(bse, q, Val(false))
		bse.Kernal(q)
	else
		kq_kindex, kΓq_kΓindex = _BSE_preprocess_eleband_qnotΓ!(bse, q, Val(true))
		bse.Kernal(kq_kindex, kΓq_kΓindex, q)
	end
	return nothing
end
function _BSE_preprocess_q!(bse::BSEgeneral, q, ::Val{false})
	q = _BSE_preprocess_eleband_q!(bse, q, Val(false))
	bse.Kernal(q)
	return nothing
end
function _BSE_preprocess_eleband_q!(bse::BSEgeneral, q, ::Val{true})
	if norm(q) < 1e-8
		q = _BSE_preprocess_eleband_q!(bse, q, Val(false))
	else
		_BSE_preprocess_eleband_qnotΓ!(bse, q, Val(true))
	end
	return q
end
function _BSE_preprocess_eleband_q!(bse::BSEgeneral, q, ::Val{false})
	# We can't calculate Γ directly.
	if norm(q) < 1e-8
		q = _BSE_shiftΓ(bse.TB.period)
	end
	bandkq = BAND(map(k -> k + q, bse.kgrid), bse.TB; vector = true)

	bse.bandkq .= bandkq
	return q
end
function _BSE_preprocess_eleband_qnotΓ!(bse::BSEgeneral, q, ::Val{true})
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

function _BSE_Hamiltonian!(bse::BSEgeneral, H)
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

			H[i, j] = Kᵈ + Kˣ
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
			H[i, i] = Δϵ + Kᵈ + Kˣ
		end
	end
	wait.(tasks)

	return Hermitian(H, :U)
end

"""
	_uijR_ψvck(bse::BSEgeneral, η)
	return an instance, which is used to extract periodic part of exciton bloch wave function.
"""
function _uijR_ψvck(bse::BSEgeneral, η)

	norb = length(bse.TB.orb_location)

	phase = _uijR_ψvck_phase(η, bse.unitcell, bse.TB.orb_location)

	vckmap = bse.vckmap
	ijRmap = bse.ijRmap

	Nvck = length(vckmap)
	NijR = length(ijRmap)

	ekR = [cis(2π * (k ⋅ R)) for R in bse.unitcell, k in bse.kgrid]
	ekqR = similar(ekR)

	UU = Array{ComplexF64}(undef, norb, norb, Nvck)
	BM = Matrix{ComplexF64}(undef, NijR, Nvck)

	return _uijR_ψvck_oneη(Nvck, NijR, vckmap, ijRmap, phase, ekR, ekqR, UU, BM)
end

function (bse::BSEgeneral)(::Val{:spinmat_vck}, q::ReducedCoordinates, spinmat_ik::AbstractMatrix{<:Number})
	norb = numorb(bse.TB)
	size(spinmat_ik) == (norb, norb) || error("Wrong spinmat_ik, its size should be (norb, norb)!")

	_BSE_preprocess_eleband_q!(bse, q, Val(bse.isqgrid))

	vckmap = bse.vckmap
	bandk = bse.bandk
	bandkq = bse.bandkq

	N = length(vckmap)
	spinmat_vck = Matrix{ComplexF64}(undef, N, N)
	tasks = Vector{Task}(undef, Int(N * (N + 1) / 2))
	n = 0
	for j in 1:N, i in 1:j
		n += 1
		tasks[n] = Threads.@spawn begin
			(v′, c′, k′) = vckmap[i]
			(v, c, k) = vckmap[j]

			#TODO The second term remain to be discussed.
			spinmat_vck[i, j] =
				bandkq[k′].vectors[:, c′] ⋅ (spinmat_ik * bandkq[k].vectors[:, c]) -
				conj(bandk[k′].vectors[:, v′]) ⋅ (spinmat_ik * conj(bandk[k].vectors[:, v]))
		end
	end
	wait.(tasks)

	return Hermitian(spinmat_vck, :U)
end
function (bse::BSEgeneral)(::Val{:spinmat_vck}, q::ReducedCoordinates, upindex::AbstractVector{<:Integer}, dnindex::AbstractVector{<:Integer} = setdiff(1:numorb(bse.TB), upindex))
	spinmat_ik = bse.TB(Val(:spinmat), upindex, dnindex)
	return bse(Val(:spinmat_vck), q, spinmat_ik)
end
function (bse::BSEgeneral)(::Val{:spinmat_ijR}, upindex, dnindex = setdiff(1:numorb(bse.TB), upindex))
	# return the spin matrix of exciton state
	return bse(H, q)
end
