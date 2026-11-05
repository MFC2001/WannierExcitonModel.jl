"""
	_uijR_phase{<: Real, <: Real}

Used to calculate the phase in excitonic u.

"""
struct _uijR_phase{qT<:Real,λT<:Real,OT<:Real}
    q::Base.RefValue{ReducedCoordinates{qT}}
    λ::Base.RefValue{λT}
    Rgrid::Vector{ReducedCoordinates{Int}}
    ndegen::Vector{Int}
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

    map(Ri -> phase.qR[Ri] = 2π * (q ⋅ phase.Rgrid[Ri]), eachindex(phase.Rgrid))
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

    map(Ri -> phase.qR[Ri] = 2π * (q ⋅ phase.Rgrid[Ri]), eachindex(phase.Rgrid))
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
    return _uijR_phase(q, λ, bse.Rgrid, bse.ndegen, bse.orb_location, bse.ijRmap)
end
function _uijR_phase(q::ReducedCoordinates{<:Real}, λ::Real,
    Rgrid::Vector{ReducedCoordinates{Int}}, ndegen::AbstractVector{<:Integer},
    orblocat::Vector{<:ReducedCoordinates{<:Real}}, ijRmap::ijRMap)

    nR = length(Rgrid)
    norb = length(orblocat)

    qR = map(x -> 2 * (q ⋅ x), Rgrid)
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

    return _uijR_phase(Ref(q), Ref(λ), Rgrid, ndegen, orblocat, ijRmap, qR, qorb, eqR, eqorb₀, phase)
end
