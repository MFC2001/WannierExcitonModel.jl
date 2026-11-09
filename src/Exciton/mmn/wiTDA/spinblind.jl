"""
    (bse::BSEspinblind_wiTDA)(:mmn, uband::AbstractVector{<:BSE_uEigen}, nnkpts::AbstractMatrix{<:Integer}; λ::Real=1//2, bandindex=nothing)

Compute the exciton overlap matrix, used for constructing maximally localized exciton Wannier functions.

`uband` must contain the exciton eigenstates on the full exciton momentum grid. 
`nnkpts` is a neighbor list specifying the neighboring k-points for each momentum point on the grid. 
`λ` is the weight of the exciton center position. 
`bandindex` selects the band indices involved in the calculation; if `nothing`, all bands are included.
"""
function (bse::BSEspinblind_wiTDA)(::Val{:mmn}, uband::AbstractVector{<:BSE_uEigen}, nnkpts::AbstractMatrix{<:Integer};
    λ::Real=1//2, bandindex=nothing)
    # uband
    if isnothing(bandindex)
        uband = map(u->u(λ).vectors, uband)
    else
        uband = map(u->u(λ, bandindex).vectors, uband)
    end
    # the center of ijR
    NijR = length(bse.ijR)
    orblocat = bse.orb_location
    center_ijR = Vector{ReducedCoordinates{Float64}}(undef, NijR)
    δ = 1 - λ
    for idx in Base.OneTo(NijR)
        (i, j, iR) = bse.ijR[idx]
        center_ijR[idx] = δ * orblocat[j] + λ * (bse.Rgrid[iR] + orblocat[i])
    end

    nband = size(uband[1], 2)
    nqb = size(nnkpts, 2)
    M = Array{ComplexF64}(undef, nband, nband, nqb)
    Threads.@threads for i in 1:nqb
        q = nnkpts[1, i]
        qb = nnkpts[2, i]
        G = ReducedCoordinates{Int}(nnkpts[3:5, i])
        uband_q = uband[q]'
        if iszero(G)
            uband_qb = uband[qb]
        else
            G_phase = map(r->cispi(-2 * (G ⋅ r)), center_ijR)
            uband_qb = G_phase .* uband[qb]
        end
        mul!(view(M, :, :, i), uband_q, uband_qb)
    end

    return M
end
