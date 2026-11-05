function (bse::BSEspinblind_wiTDA)(::Val{:amn}, )
    error("To be continued!")
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
