export BSE_uEigen
"""
	abstract type BSE_uEigen end

Used to calculate the periodic part of excitonic Bloch function.

	(u::BSE_uEigen)(λ::Real) -> Eigen

Input `λ`, return the numerical vector of periodic part.

	BSE_uEigen(bse::Union{BSESU2, BSEgeneral}, q::ReducedCoordinates{<:Real}, eigen_vck::Eigen) -> BSE_uEigen
	BSE_uEigen(bse::BSESz, q::ReducedCoordinates{<:Real}, eigen_vck::Eigen, Sz::Integer) -> BSE_uEigen

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
abstract type BSE_uEigen end
struct BSE_uEigen_spinblind{eigenT <: Eigen, qT <: Real} <: BSE_uEigen
	q::ReducedCoordinates{qT}
	eigen::eigenT
	double_qR::Vector{qT}
	double_qorb::Vector{Float64}
	ijRmap::ijRMap
end
function (u::BSE_uEigen_spinblind)(λ::Real)

	eqR = cispi.(-λ * u.double_qR)
	δ = 1 - λ
	eqorb₀ = [cispi(-λ * x - δ * y) for x in u.double_qorb, y in u.double_qorb]

	NijR = length(u.ijRmap)
	eigen_vectors = similar(u.eigen.vectors)
	for ijR in Base.OneTo(NijR)
		(i, j, R) = u.ijRmap[ijR]
		eigen_vectors[ijR, :] .= u.eigen.vectors[ijR, :] .* (eqorb₀[i, j] * eqR[R])
	end

	return Eigen(deepcopy(u.eigen.values), eigen_vectors)
end
function BSE_uEigen(bse::Union{BSESU2, BSEgeneral}, q::AbstractVector{<:Real}, eigen_vck::Eigen)

	q = ReducedCoordinates(q)

	nk = length(bse.kgrid)
	nR = length(bse.unitcell)
	e_kR = Matrix{Float64}(undef, nR, nk)
	Threads.@threads for k in Base.OneTo(nk)
		for R in Base.OneTo(nR)
			e_kR[R, k] = cispi(-2 * (bse.kgrid[k] ⋅ bse.unitcell[R]))
		end
	end
	# e_kR = [cispi(-2 * (k ⋅ R)) for R in bse.unitcell, k in bse.kgrid]
	return BSE_uEigen(bse, q, eigen_vck, e_kR)
end
function BSE_uEigen(bse::Union{BSESU2, BSEgeneral}, q::ReducedCoordinates{<:Real}, eigen_vck::Eigen, e_kR::AbstractMatrix{<:Number})

	norb = numorb(bse.TB)
	Nvck = length(bse.vckmap)
	NijR = length(bse.ijRmap)

	UU = Array{ComplexF64}(undef, norb, norb, Nvck)
	Threads.@threads for vck in Base.OneTo(Nvck)
		(v, c, k) = bse.vckmap[vck]
		for j in 1:norb, i in 1:norb
			UU[i, j, vck] = conj(bse.bandk[k].vectors[i, v]) * bse.bandkq[k].vectors[j, c]
		end
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

	double_qR = map(x -> 2 * (q ⋅ x), bse.unitcell)
	double_qorb = map(x -> 2 * (q ⋅ x), bse.TB.orb_location)

	return BSE_uEigen_spinblind(q, eigen_ijR, double_qR, double_qorb, bse.ijRmap)
end
function BSE_uEigen(bse::BSESz, q::ReducedCoordinates{<:Real}, eigen_vck::Eigen, e_kR::AbstractMatrix{<:Number}, Sz::Integer)
	return BSE_uEigen(bse, q, eigen_vck, e_kR, Val(Sz))
end
struct BSE_uEigen_Sz_0{eigenT <: Eigen, qT <: Real} <: BSE_uEigen
	q::ReducedCoordinates{qT}
	eigen::eigenT
	double_qR::Vector{qT}
	double_qorb_up::Vector{Float64}
	double_qorb_dn::Vector{Float64}
	ijRmap_uu::ijRMap
	ijRmap_dd::ijRMap
end
function (u::BSE_uEigen_Sz_0)(λ::Real)

	eqR = cispi.(-λ * u.double_qR)
	δ = 1 - λ
	eqorb₀_up = [cispi(-λ * x - δ * y) for x in u.double_qorb_up, y in u.double_qorb_up]
	eqorb₀_dn = [cispi(-λ * x - δ * y) for x in u.double_qorb_dn, y in u.double_qorb_dn]

	NijR_uu = length(u.ijRmap_uu)
	NijR_dd = length(u.ijRmap_dd)
	eigen_vectors = similar(u.eigen.vectors)
	for ijR in Base.OneTo(NijR_uu)
		(i, j, R) = u.ijRmap_uu[ijR]
		eigen_vectors[ijR, :] .= u.eigen.vectors[ijR, :] .* (eqorb₀_up[i, j] * eqR[R])
	end
	for ijR in Base.OneTo(NijR_dd)
		(i, j, R) = u.ijRmap_dd[ijR]
		ijR_uu_dd = ijR + NijR_uu
		eigen_vectors[ijR_uu_dd, :] .= u.eigen.vectors[ijR_uu_dd, :] .* (eqorb₀_dn[i, j] * eqR[R])
	end

	return Eigen(deepcopy(u.eigen.values), eigen_vectors)
end
function BSE_uEigen(bse::BSESz, q::ReducedCoordinates{<:Real}, eigen_vck::Eigen, e_kR::AbstractMatrix{<:Number}, ::Val{0})

	nk = length(bse.kgrid)
	inversesqrtnk = 1 / sqrt(nk)

	# uu
	norb_up = numorb(bse.TB_up)
	Nvck_uu = length(bse.vckmap_uu)
	UU_uu = Array{ComplexF64}(undef, norb_up, norb_up, Nvck_uu)
	Threads.@threads for vck in Base.OneTo(Nvck_uu)
		(v, c, k) = bse.vckmap_uu[vck]
		for j in 1:norb_up, i in 1:norb_up
			UU_uu[i, j, vck] = conj(bse.bandk_up[k].vectors[i, v]) * bse.bandkq_up[k].vectors[j, c]
		end
	end
	NijR_uu = length(bse.ijRmap_uu)
	BM_uu = Matrix{ComplexF64}(undef, NijR_uu, Nvck_uu)
	Threads.@threads for vck in Base.OneTo(Nvck_uu)
		(_, _, k) = bse.vckmap_uu[vck]
		for ijR in Base.OneTo(NijR_uu)
			(i, j, R) = bse.ijRmap_uu[ijR]
			BM_uu[ijR, vck] = inversesqrtnk * UU_uu[i, j, vck] * e_kR[R, k]
		end
	end

	# dd
	norb_dn = numorb(bse.TB_dn)
	Nvck_dd = length(bse.vckmap_dd)
	UU_dd = Array{ComplexF64}(undef, norb_dn, norb_dn, Nvck_dd)
	Threads.@threads for vck in Base.OneTo(Nvck_dd)
		(v, c, k) = bse.vckmap_dd[vck]
		for j in 1:norb_dn, i in 1:norb_dn
			UU_dd[i, j, vck] = conj(bse.bandk_dn[k].vectors[i, v]) * bse.bandkq_dn[k].vectors[j, c]
		end
	end
	NijR_dd = length(bse.ijRmap_dd)
	BM_dd = Matrix{ComplexF64}(undef, NijR_dd, Nvck_dd)
	Threads.@threads for vck in Base.OneTo(Nvck_dd)
		(_, _, k) = bse.vckmap_dd[vck]
		for ijR in Base.OneTo(NijR_dd)
			(i, j, R) = bse.ijRmap_dd[ijR]
			BM_dd[ijR, vck] = inversesqrtnk * UU_dd[i, j, vck] * e_kR[R, k]
		end
	end

	eigen_ijR_vectors = Matrix{ComplexF64}(undef, NijR_uu + NijR_dd, size(eigen_vck.vectors, 2))
	mul!(view(eigen_ijR_vectors, 1:NijR_uu, :), BM_uu, view(eigen_vck.vectors, 1:Nvck_uu, :))
	mul!(view(eigen_ijR_vectors, NijR_uu+1:NijR_uu+NijR_dd, :), BM_dd, view(eigen_vck.vectors, Nvck_uu+1:Nvck_uu+Nvck_dd, :))

	eigen_ijR = Eigen(deepcopy(eigen_vck.values), eigen_ijR_vectors)

	double_qR = map(x -> 2 * (q ⋅ x), bse.unitcell)
	double_qorb_up = map(x -> 2 * (q ⋅ x), bse.TB_up.orb_location)
	double_qorb_dn = map(x -> 2 * (q ⋅ x), bse.TB_dn.orb_location)

	return BSE_uEigen_Sz_0(q, eigen_ijR, double_qR, double_qorb_up, double_qorb_dn, bse.ijRmap_uu, bse.ijRmap_dd)
end
struct BSE_uEigen_Sz_pn1{eigenT <: Eigen, qT <: Real} <: BSE_uEigen
	q::ReducedCoordinates{qT}
	eigen::eigenT
	double_qR::Vector{qT}
	double_qorb_v::Vector{Float64}
	double_qorb_c::Vector{Float64}
	ijRmap::ijRMap
end
function (u::BSE_uEigen_Sz_pn1)(λ::Real)

	eqR = cispi.(-λ * u.double_qR)
	δ = 1 - λ
	eqorb₀_up = [cispi(-λ * x - δ * y) for x in u.double_qorb_v, y in u.double_qorb_c]

	NijR = length(u.ijRmap)
	eigen_vectors = similar(u.eigen.vectors)
	for ijR in Base.OneTo(NijR)
		(i, j, R) = u.ijRmap[ijR]
		eigen_vectors[ijR, :] .= u.eigen.vectors[ijR, :] .* (eqorb₀_up[i, j] * eqR[R])
	end

	return Eigen(deepcopy(u.eigen.values), eigen_vectors)
end
function BSE_uEigen(bse::BSESz, q::ReducedCoordinates{<:Real}, eigen_vck::Eigen, e_kR::AbstractMatrix{<:Number}, ::Val{1})

	nk = length(bse.kgrid)
	inversesqrtnk = 1 / sqrt(nk)

	# du
	norb_up = numorb(bse.TB_up)
	norb_dn = numorb(bse.TB_dn)
	Nvck_du = length(bse.vckmap_du)
	UU_du = Array{ComplexF64}(undef, norb_dn, norb_up, Nvck_du)
	Threads.@threads for vck in Base.OneTo(Nvck_du)
		(v, c, k) = bse.vckmap_du[vck]
		for j in 1:norb_up, i in 1:norb_dn
			UU_du[i, j, vck] = conj(bse.bandk_dn[k].vectors[i, v]) * bse.bandkq_up[k].vectors[j, c]
		end
	end
	NijR_du = length(bse.ijRmap_du)
	BM_du = Matrix{ComplexF64}(undef, NijR_du, Nvck_du)
	Threads.@threads for vck in Base.OneTo(Nvck_du)
		(_, _, k) = bse.vckmap_du[vck]
		for ijR in Base.OneTo(NijR_du)
			(i, j, R) = bse.ijRmap_du[ijR]
			BM_du[ijR, vck] = inversesqrtnk * UU_du[i, j, vck] * e_kR[R, k]
		end
	end

	eigen_ijR_vectors = BM_du * eigen_vck.vectors
	eigen_ijR = Eigen(deepcopy(eigen_vck.values), eigen_ijR_vectors)

	double_qR = map(x -> 2 * (q ⋅ x), bse.unitcell)
	double_qorb_up = map(x -> 2 * (q ⋅ x), bse.TB_up.orb_location)
	double_qorb_dn = map(x -> 2 * (q ⋅ x), bse.TB_dn.orb_location)

	return BSE_uEigen_Sz_pn1(q, eigen_ijR, double_qR, double_qorb_dn, double_qorb_up, bse.ijRmap_du)
end
function BSE_uEigen(bse::BSESz, q::ReducedCoordinates{<:Real}, eigen_vck::Eigen, e_kR::AbstractMatrix{<:Number}, ::Val{-1})

	nk = length(bse.kgrid)
	inversesqrtnk = 1 / sqrt(nk)

	# ud
	norb_up = numorb(bse.TB_up)
	norb_dn = numorb(bse.TB_dn)
	Nvck_ud = length(bse.vckmap_ud)
	UU_ud = Array{ComplexF64}(undef, norb_up, norb_dn, Nvck_ud)
	Threads.@threads for vck in Base.OneTo(Nvck_ud)
		(v, c, k) = bse.vckmap_ud[vck]
		for j in 1:norb_dn, i in 1:norb_up
			UU_ud[i, j, vck] = conj(bse.bandk_up[k].vectors[i, v]) * bse.bandkq_dn[k].vectors[j, c]
		end
	end
	NijR_ud = length(bse.ijRmap_ud)
	BM_ud = Matrix{ComplexF64}(undef, NijR_ud, Nvck_ud)
	Threads.@threads for vck in Base.OneTo(Nvck_ud)
		(_, _, k) = bse.vckmap_ud[vck]
		for ijR in Base.OneTo(NijR_ud)
			(i, j, R) = bse.ijRmap_ud[ijR]
			BM_ud[ijR, vck] = inversesqrtnk * UU_ud[i, j, vck] * e_kR[R, k]
		end
	end

	eigen_ijR_vectors = BM_ud * eigen_vck.vectors
	eigen_ijR = Eigen(deepcopy(eigen_vck.values), eigen_ijR_vectors)

	double_qR = map(x -> 2 * (q ⋅ x), bse.unitcell)
	double_qorb_up = map(x -> 2 * (q ⋅ x), bse.TB_up.orb_location)
	double_qorb_dn = map(x -> 2 * (q ⋅ x), bse.TB_dn.orb_location)

	return BSE_uEigen_Sz_pn1(q, eigen_ijR, double_qR, double_qorb_up, double_qorb_dn, bse.ijRmap_ud)
end
"""
	_uijR_phase{<: Real, <: Real}

Used to calculate the phase in excitonic u.

"""
struct _uijR_phase{qT <: Real, λT <: Real, OT <: Real}
	q::Base.RefValue{ReducedCoordinates{qT}}
	λ::Base.RefValue{λT}
	unitcell::Vector{ReducedCoordinates{Int}}
	orblocat::Vector{ReducedCoordinates{OT}}
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
		phase.eqR[Ri] = cispi(-λ * qR)
	end
	for j in Base.OneTo(norb), i in Base.OneTo(norb)
		phase.eqorb₀[i, j] = cispi(-δ * phase.qorb[i] - λ * phase.qorb[j])
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

	map(Ri -> phase.qR[Ri] = 2π * (q ⋅ phase.unitcell[Ri]), eachindex(phase.unitcell))
	map(i -> phase.qorb[i] = 2π * (q ⋅ phase.orblocat[i]), eachindex(phase.orblocat))

	for (Ri, qR) in enumerate(phase.qR)
		phase.eqR[Ri] = cispi(-phase.λ[] * qR)
	end
	δ = 1 - phase.λ[]
	for j in Base.OneTo(norb), i in Base.OneTo(norb)
		phase.eqorb₀[i, j] = cispi(-δ * phase.qorb[i] - phase.λ[] * phase.qorb[j])
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

	map(Ri -> phase.qR[Ri] = 2π * (q ⋅ phase.unitcell[Ri]), eachindex(phase.unitcell))
	map(i -> phase.qorb[i] = 2π * (q ⋅ phase.orblocat[i]), eachindex(phase.orblocat))

	for (Ri, qR) in enumerate(phase.qR)
		phase.eqR[Ri] = cispi(-λ * qR)
	end
	for j in Base.OneTo(norb), i in Base.OneTo(norb)
		phase.eqorb₀[i, j] = cispi(-δ * phase.qorb[i] - λ * phase.qorb[j])
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
	unitcell::Vector{ReducedCoordinates{Int}}, orblocat::Vector{<:ReducedCoordinates{<:Real}}, ijRmap::ijRMap)

	nR = length(unitcell)
	norb = length(orblocat)

	qR = map(x -> 2 * (q ⋅ x), unitcell)
	qorb = map(x -> 2 * (q ⋅ x), orblocat)

	eqR = cispi.(-λ * qR)

	δ = 1 - λ
	eqorb₀ = Matrix{ComplexF64}(undef, norb, norb)
	for j in Base.OneTo(norb), i in Base.OneTo(norb)
		eqorb₀[i, j] = cispi(-δ * qorb[i] - λ * qorb[j])
	end

	NijR = length(ijRmap)
	phase = Vector{ComplexF64}(undef, NijR)
	for ijR in Base.OneTo(NijR)
		(i, j, R) = ijRmap[ijR]
		phase[ijR] = eqorb₀[i, j] * eqR[R]
	end

	return _uijR_phase(Ref(q), Ref(λ), unitcell, orblocat, ijRmap, qR, qorb, eqR, eqorb₀, phase)
end
