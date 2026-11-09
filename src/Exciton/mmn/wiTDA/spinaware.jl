"""
    (bse::BSEspinaware_wiTDA)(:mmn, uband::AbstractVector{<:BSE_uEigen}, nnkpts::AbstractMatrix{<:Integer}, Sz::Integer; λ::Real=1//2, bandindex=nothing)

Compute the exciton overlap matrix, used for constructing maximally localized exciton Wannier functions.

`uband` must contain the exciton eigenstates on the full exciton momentum grid, and `Sz` is the z-component of the exciton's total spin.
`nnkpts` is a neighbor list specifying the neighboring k-points for each momentum point on the grid. 
`λ` is the weight of the exciton center position. 
`bandindex` selects the band indices involved in the calculation; if `nothing`, all bands are included.
"""
function (bse::BSEspinaware_wiTDA)(::Val{:mmn}, uband::AbstractVector{<:BSE_uEigen}, nnkpts::AbstractMatrix{<:Integer}, Sz::Integer;
    λ::Real=1//2, bandindex=nothing)
    # uband
    if isnothing(bandindex)
        uband = map(u->u(λ).vectors, uband)
    else
        uband = map(u->u(λ, bandindex).vectors, uband)
    end
    center_ijR = _mmn_center_BSEspinaware_wiTDA(bse, λ, Val(Sz))

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
function _mmn_center_BSEspinaware_wiTDA(bse::BSEspinaware_wiTDA, λ::Real, ::Val{0})
    # the center of ijR
    NijR_uu = length(bse.ijR_uu)
    NijR_dd = length(bse.ijR_dd)
    orblocat_up = bse.orb_location_up
    orblocat_dn = bse.orb_location_dn
    center_ijR = Vector{ReducedCoordinates{Float64}}(undef, NijR_uu + NijR_dd)
    δ = 1 - λ
    for idx in Base.OneTo(NijR_uu)
        (i, j, iR) = bse.ijR_uu[idx]
        center_ijR[idx] = δ * orblocat_up[j] + λ * (bse.Rgrid[iR] + orblocat_up[i])
    end
    for idx in Base.OneTo(NijR_dd)
        (i, j, iR) = bse.ijR_dd[idx]
        center_ijR[idx+NijR_uu] = δ * orblocat_dn[j] + λ * (bse.Rgrid[iR] + orblocat_dn[i])
    end
    return center_ijR
end
function _mmn_center_BSEspinaware_wiTDA(bse::BSEspinaware_wiTDA, λ::Real, ::Val{-1})
    # the center of ijR
    NijR_ud = length(bse.ijR_ud)
    orblocat_up = bse.orb_location_up
    orblocat_dn = bse.orb_location_dn
    center_ijR = Vector{ReducedCoordinates{Float64}}(undef, NijR_ud)
    δ = 1 - λ
    for idx in Base.OneTo(NijR_ud)
        (i, j, iR) = bse.ijR_ud[idx]
        center_ijR[idx] = δ * orblocat_dn[j] + λ * (bse.Rgrid[iR] + orblocat_up[i])
    end
    return center_ijR
end
function _mmn_center_BSEspinaware_wiTDA(bse::BSEspinaware_wiTDA, λ::Real, ::Val{1})
    # the center of ijR
    NijR_du = length(bse.ijR_du)
    orblocat_up = bse.orb_location_up
    orblocat_dn = bse.orb_location_dn
    center_ijR = Vector{ReducedCoordinates{Float64}}(undef, NijR_du)
    δ = 1 - λ
    for idx in Base.OneTo(NijR_du)
        (i, j, iR) = bse.ijR_du[idx]
        center_ijR[idx] = δ * orblocat_up[j] + λ * (bse.Rgrid[iR] + orblocat_dn[i])
    end
    return center_ijR
end
