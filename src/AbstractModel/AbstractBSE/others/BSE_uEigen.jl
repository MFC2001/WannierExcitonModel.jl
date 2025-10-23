"""
	BSE_uEigen{<: Eigen}

Used to calculate the periodic part of excitonic Bloch function.

	(u::BSE_uEigen)(λ::Real) -> Eigen

Input `λ`, return the numerical vector of periodic part.

	BSE_uEigen(bse::AbstractBSE, q::ReducedCoordinates{<:Real}, eigen_vck::Eigen) -> BSE_uEigen

Create an object of `BSE_uEigen`.
`eigen_vck` contains the excitonic states.
!!!note

	This method will use the electronic Bloch functions, which are stored in `AbstractBSE` after creating a BSE Hamiltonian.
	Ensuring that the last progress of `bse` is about `q`.

```julia
julia> H = bse(q)
julia> eigen_vck = eigen(H)
julia> u = BSE_uEigen(bse, q, eigen_vck)
julia> u_λ = u(λ)
```
"""
struct BSE_uEigen{eigenT <: Eigen}
	eigen::eigenT
	qR::Vector{Float64}
	qorb::Vector{Float64}
	ijRmap::ijRMap
end
function (u::BSE_uEigen)(λ::Real)

	eqR = cis.(-λ * u.qR)
	norb = length(u.qorb)
	eqorb₀ = [cis(-δ * u.qorb[i] - λ * u.qorb[j]) for i in Base.OneTo(norb), j in Base.OneTo(norb)]

	NijR = length(u.ijRmap)
	eigen_vectors = deepcopy(u.eigen.vectors)
	for ijR in Base.OneTo(NijR)
		(i, j, R) = u.ijRmap[ijR]
		eigen_vectors[ijR, :] .*= eqorb₀[i, j] * eqR[R]
	end

	return Eigen(deepcopy(u.eigen.values), eigen_vectors)
end
function BSE_uEigen(bse::AbstractBSE, q::ReducedCoordinates{<:Real}, eigen_vck::Eigen)

	norb = numorb(bse.TB)
	Nvck = length(bse.vckmap)
	NijR = length(ijRmap)

	UU = Array{ComplexF64}(undef, norb, norb, Nvck)
	Threads.@threads for vck in Base.OneTo(Nvck)
		(v, c, k) = bse.vckmap[vck]
		UU[:, :, vck] .= bse.bandkq[k].vectors[:, c] * (bse.bandk[k].vectors[:, v])'
	end

	nk = length(bse.kgrid)
	nR = length(bse.unitcell)
	e_kR = Matrix{Float64}(undef, nR, nk)
	Threads.@threads for k in Base.OneTo(nk)
		for R in Base.OneTo(nR)
			e_kR[R, k] = cis(-2π * (bse.kgrid[k] ⋅ bse.unitcell[R]))
		end
	end
	# e_kR = [cis(-2π * (k ⋅ R)) for R in bse.unitcell, k in bse.kgrid]


	inversesqrtnk = 1 / sqrt(nk)

	BM = Matrix{ComplexF64}(undef, NijR, Nvck)
	Threads.@threads for vck in Base.OneTo(Nvck)
		(_, _, k) = bse.vckmap[vck]
		for ijR in Base.OneTo(NijR)
			(i, j, R) = bse.ijRmap[ijR]
			BM[ijR, vck] = inversesqrtnk * UU[i, j, vck] * e_kR[R, k]
		end
	end

	eigen_ijR_vectors = BM * eigen_vck.vectors
	eigen_ijR = Eigen(deepcopy(eigen_vck.values), eigen_ijR_vectors)

	qR = map(x -> 2π * (q ⋅ x), bse.unitcell)
	qorb = map(x -> 2π * (q ⋅ x), bse.TB.orb_location)

	return BSE_uEigen(eigen_ijR, qR, qorb, bse.ijRmap)
end
function BSE_uEigen(bse::AbstractBSE, q::ReducedCoordinates{<:Real}, eigen_vck::Eigen, e_kR::Matrix{<:Number})

	norb = numorb(bse.TB)
	Nvck = length(bse.vckmap)
	NijR = length(ijRmap)

	UU = Array{ComplexF64}(undef, norb, norb, Nvck)
	Threads.@threads for vck in Base.OneTo(Nvck)
		(v, c, k) = bse.vckmap[vck]
		UU[:, :, vck] .= bse.bandkq[k].vectors[:, c] * (bse.bandk[k].vectors[:, v])'
	end

	nk = length(bse.kgrid)
	inversesqrtnk = 1 / sqrt(nk)

	BM = Matrix{ComplexF64}(undef, NijR, Nvck)
	Threads.@threads for vck in Base.OneTo(Nvck)
		(_, _, k) = bse.vckmap[vck]
		for ijR in Base.OneTo(NijR)
			(i, j, R) = bse.ijRmap[ijR]
			BM[ijR, vck] = inversesqrtnk * UU[i, j, vck] * e_kR[R, k]
		end
	end

	eigen_ijR_vectors = BM * eigen_vck.vectors
	eigen_ijR = Eigen(deepcopy(eigen_vck.values), eigen_ijR_vectors)

	qR = map(x -> 2π * (q ⋅ x), bse.unitcell)
	qorb = map(x -> 2π * (q ⋅ x), bse.TB.orb_location)

	return BSE_uEigen(eigen_ijR, qR, qorb, bse.ijRmap)
end
"""
	_uijR_phase{<: Real, <: Real}

Used to calculate the phase in excitonic u.

"""
struct _uijR_phase
	q::Ref{ReducedCoordinates}
	λ::Ref{Real}
	R::Vector{ReducedCoordinates{Int}}
	orblocat::Vector{ReducedCoordinates}
	ijRmap::ijRMap
	qR::Vector{Float64}
	qorb::Vector{Float64}
	eqR::Vector{ComplexF64}
	eqorb₀::Matrix{ComplexF64}
	phase::Vector{ComplexF64}
end
function (phase::_uijR_phase)()
	return phase.phase
end
function (phase::_uijR_phase)(λ::Real)

	phase.λ[] = λ

	norb = length(phase.orb)
	δ = 1 - λ

	# for j in Base.OneTo(norb), i in Base.OneTo(norb)
	# 	phase.orb₀[i, j] = δ * phase.orblocat[i] + η * phase.orblocat[j]
	# end

	for (Ri, qR) in enumerate(phase.qR)
		phase.eqR[Ri] = cis(-λ * qR)
	end
	for j in Base.OneTo(norb), i in Base.OneTo(norb)
		phase.eqorb₀[i, j] = cis(-δ * phase.qorb[i] - λ * phase.qorb[j])
	end

	NijR = length(phase.ijRmap)
	for ijR in Base.OneTo(NijR)
		(i, j, R) = phase.ijRmap[ijR]
		phase.phase[ijR] = phase.eqorb₀[i, j] * phase.eqR[R]
	end

	return phase.phase
end
function (phase::_uijR_phase)(q::ReducedCoordinates{<:Real})
	phase.q[] = q

	map(Ri -> phase.qR[Ri] = 2π * (q ⋅ phase.R[Ri]), eachindex(phase.R))
	map(i -> phase.qorb[i] = 2π * (q ⋅ phase.orblocat[i]), eachindex(phase.orblocat))

	for (Ri, qR) in enumerate(phase.qR)
		phase.eqR[Ri] = cis(-phase.λ[] * qR)
	end
	δ = 1 - phase.λ[]
	for j in Base.OneTo(norb), i in Base.OneTo(norb)
		phase.eqorb₀[i, j] = cis(-δ * phase.qorb[i] - phase.λ[] * phase.qorb[j])
	end

	NijR = length(phase.ijRmap)
	for ijR in Base.OneTo(NijR)
		(i, j, R) = phase.ijRmap[ijR]
		phase.phase[ijR] = phase.eqorb₀[i, j] * phase.eqR[R]
	end

	return phase.phase
end
function (phase::_uijR_phase)(q::ReducedCoordinates{<:Real}, λ::Real)

	phase.q[] = q
	phase.λ[] = λ

	norb = length(phase.orb)
	δ = 1 - λ

	map(Ri -> phase.qR[Ri] = 2π * (q ⋅ phase.R[Ri]), eachindex(phase.R))
	map(i -> phase.qorb[i] = 2π * (q ⋅ phase.orblocat[i]), eachindex(phase.orblocat))

	for (Ri, qR) in enumerate(phase.qR)
		phase.eqR[Ri] = cis(-λ * qR)
	end
	for j in Base.OneTo(norb), i in Base.OneTo(norb)
		phase.eqorb₀[i, j] = cis(-δ * phase.qorb[i] - λ * phase.qorb[j])
	end

	NijR = length(phase.ijRmap)
	for ijR in Base.OneTo(NijR)
		(i, j, R) = phase.ijRmap[ijR]
		phase.phase[ijR] = phase.eqorb₀[i, j] * phase.eqR[R]
	end

	return phase.phase
end
function _uijR_phase(bse::AbstractBSE, q::ReducedCoordinates{<:Real}, λ::Real)
	return _uijR_phase(q, λ, bse.unitcell, bse.TB.orb_location, bse.ijRmap)
end
function _uijR_phase(q::ReducedCoordinates{<:Real}, λ::Real,
	R::Vector{ReducedCoordinates{Int}}, orblocat::Vector{<:ReducedCoordinates{<:Real}}, ijRmap::ijRMap)

	nR = length(R)
	norb = length(orblocat)

	qR = map(x -> 2π * (q ⋅ x), R)
	qorb = map(x -> 2π * (q ⋅ x), orblocat)

	eqR = cis.(-λ * qR)

	δ = 1 - λ
	eqorb₀ = Matrix{ComplexF64}(undef, norb, norb)
	for j in Base.OneTo(norb), i in Base.OneTo(norb)
		eqorb₀[i, j] = cis(-δ * qorb[i] - λ * qorb[j])
	end

	NijR = length(ijRmap)
	phase = Vector{ComplexF64}(undef, NijR)
	for ijR in Base.OneTo(NijR)
		(i, j, R) = ijRmap[ijR]
		phase[ijR] = eqorb₀[i, j] * eqR[R]
	end

	return _uijR_phase(Ref(q), Ref(λ), R, orblocat, ijRmap, qR, qorb, eqR, eqorb₀, phase)
end
