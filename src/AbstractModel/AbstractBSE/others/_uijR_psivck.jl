"""
	_uijR_ψvck_phase{<: Real, <: Real}

Used to calculate `eqR` and `eqorb₀` for each `q`.
"""
struct _uijR_ψvck_phase{T <: Real, S <: Real}
	η::T
	invsqrtnR::Float64
	R::Vector{ReducedCoordinates{Int}}
	orb₀::Matrix{ReducedCoordinates{S}}
	eqR::Vector{ComplexF64}
	eqorb₀::Matrix{ComplexF64}
end
function (phase::_uijR_ψvck_phase)(q)
	for (Ri, R) in enumerate(phase.R)
		phase.eqR[Ri] = phase.invsqrtnR * cis(2π * phase.η * (q ⋅ R))
	end
	for I in CartesianIndices(phase.orb₀)
		phase.eqorb₀[I] = cis(-2π * (q ⋅ phase.orb₀[I]))
	end
end
function _uijR_ψvck_phase(η, R, orblocat)

	norb = length(orblocat)
	nR = length(R)

	invsqrtnR = 1 / sqrt(nR)

	δ = 1 - η
	orb₀ = map(CartesianIndices((norb, norb))) do I
		(i, j) = Tuple(I)
		(δ * orblocat[i] + η * orblocat[j])
	end

	eqorb₀ = Matrix{ComplexF64}(undef, norb, norb)
	eqR = Vector{ComplexF64}(undef, nR)

	return _uijR_ψvck_phase(η, invsqrtnR, R, orb₀, eqR, eqorb₀)
end
"""
	_uijR_ψvck_oneη

Used to calculate `BM` for each `q`.
u = BM * ψ
This is used when ηt = ηs, so we only need to calculate `BM` once.
"""
struct _uijR_ψvck_oneη{T, S}
	Nvck::Int
	NijR::Int
	vckmap::vckMap
	ijRmap::ijRMap
	phase::_uijR_ψvck_phase{T, S}
	ekR::Matrix{ComplexF64}
	ekqR::Matrix{ComplexF64}
	UU::Array{ComplexF64}
	BM::Matrix{ComplexF64}
end
function (ijRvck::_uijR_ψvck_oneη)(bandk, bandkq, q)

	Threads.@threads for vck in Base.OneTo(ijRvck.Nvck)
		(v, c, k) = ijRvck.vckmap[vck]
		ijRvck.UU[:, :, vck] .= bandkq[k].vectors[:, c] * transpose(conj(bandk[k].vectors[:, v]))
	end

	ijRvck.phase(q)

	ijRvck.ekqR .= ijRvck.phase.eqR .* ijRvck.ekR
	Threads.@threads for vck in Base.OneTo(ijRvck.Nvck)
		(_, _, k) = ijRvck.vckmap[vck]
		for ijR in Base.OneTo(ijRvck.NijR)
			(i, j, R) = ijRvck.ijRmap[ijR]
			ijRvck.BM[ijR, vck] = ijRvck.UU[i, j, vck] * ijRvck.ekqR[R, k] * ijRvck.phase.eqorb₀[i, j]
		end
	end

	return ijRvck.BM, ijRvck.BM
end
"""
	_uijR_ψvck_twoη

Used to calculate `BM` for each `q`.
u = BM * ψ
This will calculate `BMt` and `BMs` for triplet and singlet at the same time.
"""
struct _uijR_ψvck_twoη{Tt, Ts, St, Ss}
	Nvck::Int
	NijR::Int
	vckmap::vckMap
	ijRmap::ijRMap
	phase_t::_uijR_ψvck_phase{Tt, St}
	phase_s::_uijR_ψvck_phase{Ts, Ss}
	ekR::Matrix{ComplexF64}
	ekqR::Matrix{ComplexF64}
	UU::Array{ComplexF64}
	BMt::Matrix{ComplexF64}
	BMs::Matrix{ComplexF64}
end
function (ijRvck::_uijR_ψvck_twoη)(bandk, bandkq, q)

	Threads.@threads for vck in Base.OneTo(ijRvck.Nvck)
		(v, c, k) = ijRvck.vckmap[vck]
		ijRvck.UU[:, :, vck] .= bandkq[k].vectors[:, c] * transpose(conj(bandk[k].vectors[:, v]))
	end

	ijRvck.phase_t(q)

	ijRvck.ekqR .= ijRvck.phase_t.eqR .* ijRvck.ekR
	Threads.@threads for vck in Base.OneTo(ijRvck.Nvck)
		(_, _, k) = ijRvck.vckmap[vck]
		for ijR in Base.OneTo(ijRvck.NijR)
			(i, j, R) = ijRvck.ijRmap[ijR]
			ijRvck.BMt[ijR, vck] = ijRvck.UU[i, j, vck] * ijRvck.ekqR[R, k] * ijRvck.phase_t.eqorb₀[i, j]
		end
	end

	ijRvck.phase_s(q)

	ijRvck.ekqR .= ijRvck.phase_s.eqR .* ijRvck.ekR
	Threads.@threads for vck in Base.OneTo(ijRvck.Nvck)
		(_, _, k) = ijRvck.vckmap[vck]
		for ijR in Base.OneTo(ijRvck.NijR)
			(i, j, R) = ijRvck.ijRmap[ijR]
			ijRvck.BMs[ijR, vck] = ijRvck.UU[i, j, vck] * ijRvck.ekqR[R, k] * ijRvck.phase_s.eqorb₀[i, j]
		end
	end

	return ijRvck.BMt, ijRvck.BMs
end

# """
# 	_uijR_ψvck_spinful
# """
# struct _uijR_ψvck_spinful{T, S}
# 	Nvck::Int
# 	NijR::Int
# 	vckmap::vckMap
# 	ijRmap::ijRMap
# 	phase::_uijR_ψvck_phase{T, S}
# 	ekR::Matrix{ComplexF64}
# 	ekqR::Matrix{ComplexF64}
# 	UU::Array{ComplexF64}
# 	BM::Matrix{ComplexF64}
# end
# function (ijRvck::_uijR_ψvck_spinful)(bandk, bandkq, q)

# 	Threads.@threads for vck in Base.OneTo(ijRvck.Nvck)
# 		(v, c, k) = ijRvck.vckmap[vck]
# 		ijRvck.UU[:, :, vck] .= bandkq[k].vectors[:, c] * transpose(conj(bandk[k].vectors[:, v]))
# 	end

# 	ijRvck.phase(q)

# 	ijRvck.ekqR .= ijRvck.phase.eqR .* ijRvck.ekR
# 	Threads.@threads for vck in Base.OneTo(ijRvck.Nvck)
# 		(_, _, k) = ijRvck.vckmap[vck]
# 		for ijR in Base.OneTo(ijRvck.NijR)
# 			(i, j, R) = ijRvck.ijRmap[ijR]
# 			ijRvck.BM[ijR, vck] = ijRvck.UU[i, j, vck] * ijRvck.ekqR[R, k] * ijRvck.phase.eqorb₀[i, j]
# 		end
# 	end

# 	return ijRvck.BM
# end

