"""
	BSE_uEigen(bse::BSEspinaware_wiTDA, q::ReducedCoordinates{<:Real}, eigen_vck::Eigen, Sz::Integer) -> BSE_uEigen

Create an object of `BSE_uEigen`.
`eigen_vck` contains the excitonic states.
`Sz` is the excitonic Sz.

```julia
julia> H = bse(0, q)
julia> eigen_vck = eigen(H)
julia> u = bse(:periodic, q, eigen_vck, 0)
```
"""
function (bse::BSEspinaware_wiTDA)(::Val{:periodic}, q::AbstractVector{<:Real}, eigen_vck::Eigen, Sz::Integer)
    if Sz == 0
        Nvck_Sz0 = length(bse.vck_uu) + length(bse.vck_dd)
        @assert !iszero(Nvck_Sz0) "Cannot work for nonexistent Sz-space"
        size(eigen_vck.vectors, 1) == Nvck_Sz0 || throw(DimensionMismatch(
            "vector length should be Nvck=$(Nvck_Sz0), got $(size(eigen_vck.vectors, 1))"))
    elseif Sz == -1
        Nvck_Szn1 = length(bse.vck_ud)
        @assert !iszero(Nvck_Szn1) "Cannot work for nonexistent Sz-space"
        size(eigen_vck.vectors, 1) == Nvck_Szn1 || throw(DimensionMismatch(
            "vector length should be Nvck=$(length(bse.vck_ud)), got $(size(eigen_vck.vectors, 1))"))
    elseif Sz == 1
        Nvck_Szp1 = length(bse.vck_du)
        @assert !iszero(Nvck_Szp1) "Cannot work for nonexistent Sz-space"
        size(eigen_vck.vectors, 1) == length(bse.vck_du) || throw(DimensionMismatch(
            "vector length should be Nvck=$(length(bse.vck_du)), got $(size(eigen_vck.vectors, 1))"))
    else
        error("Sz must be 0, -1, or 1, got $Sz")
    end
    # e_kR
    nk = length(bse.kgrid)
    nR = length(bse.Rgrid)
    e_kR = Matrix{ComplexF64}(undef, nR, nk)
    Threads.@threads for ik in Base.OneTo(nk)
        k = bse.kgrid[ik]
        for (iR, R) in enumerate(bse.Rgrid)
            e_kR[iR, ik] = cispi(-2 * (k ⋅ R))
        end
    end
    # bandkq
    q = ReducedCoordinates(q)
    kq = map(k -> k + q, bse.kgrid)
    if Sz == 0
        if iszero(length(bse.vck_uu))
            bandkq_up = nothing
        else
            bandkq_up = BAND(kq, bse.eleH_up; vector=true)
            phase_normalize!(:positive_maxabs, bandkq_up)
        end
        if iszero(length(bse.vck_dd))
            bandkq_dn = nothing
        else
            bandkq_dn = BAND(kq, bse.eleH_dn; vector=true)
            phase_normalize!(:positive_maxabs, bandkq_dn)
        end
        return _BSE_uEigen(bse, q, eigen_vck, bandkq_up, bandkq_dn, e_kR, Val(0))
    elseif Sz == -1
        bandkq_up = nothing
        bandkq_dn = BAND(kq, bse.eleH_dn; vector=true)
        phase_normalize!(:positive_maxabs, bandkq_dn)
        return _BSE_uEigen(bse, q, eigen_vck, bandkq_dn, e_kR, Val(-1))
    else # Sz == 1
        bandkq_dn = nothing
        bandkq_up = BAND(kq, bse.eleH_up; vector=true)
        phase_normalize!(:positive_maxabs, bandkq_up)
        return _BSE_uEigen(bse, q, eigen_vck, bandkq_up, e_kR, Val(1))
    end
end
# Sz0
struct BSEspinaware_wiTDA_uEigen_Sz0{eigenT<:Eigen,qT<:Real} <: BSE_uEigen
    q::ReducedCoordinates{qT}
    eigen::eigenT
    inv_sqrt_ndegen::Vector{Float64}
    double_qR::Vector{qT}
    double_qorb_up::Vector{Float64}
    double_qorb_dn::Vector{Float64}
    ijR_uu::TripleList{Int,Int,Int}
    ijR_dd::TripleList{Int,Int,Int}
end
function (u::BSEspinaware_wiTDA_uEigen_Sz0)(λ::Real, bandindex::AbstractVector{<:Integer}=eachindex(u.eigen.values))

    eqR = cispi.(-λ * u.double_qR) .* u.inv_sqrt_ndegen
    δ = 1 - λ
    eqorb₀_up = [cispi(-λ * x - δ * y) for x in u.double_qorb_up, y in u.double_qorb_up]
    eqorb₀_dn = [cispi(-λ * x - δ * y) for x in u.double_qorb_dn, y in u.double_qorb_dn]

    NijR_uu = length(u.ijR_uu)
    NijR_dd = length(u.ijR_dd)
    eigen_vectors = Matrix{eltype(u.eigen.vectors)}(undef, NijR_uu + NijR_dd, length(bandindex))
    Threads.@threads for idx in Base.OneTo(NijR_uu)
        (i, j, iR) = u.ijR_uu[idx]
        eigen_vectors[idx, :] .= u.eigen.vectors[idx, bandindex] .* (eqorb₀_up[i, j] * eqR[iR])
    end
    Threads.@threads for idx in Base.OneTo(NijR_dd)
        (i, j, iR) = u.ijR_dd[idx]
        idx_dd = idx + NijR_uu
        eigen_vectors[idx_dd, :] .= u.eigen.vectors[idx_dd, bandindex] .* (eqorb₀_dn[i, j] * eqR[iR])
    end

    return Eigen(u.eigen.values[bandindex], eigen_vectors)
end
function _BSE_uEigen(bse::BSEspinaware_wiTDA, q::ReducedCoordinates{qT}, eigen_vck::Eigen,
    bandkq_up, bandkq_dn, e_kR, ::Val{0}) where {qT}

    BM_uu, BM_dd = _vck2ijR(bse, bandkq_up, bandkq_dn, e_kR, Val(0))

    Nvck_uu = length(bse.vck_uu)
    Nvck_dd = length(bse.vck_dd)
    vectors_ijR_uu = BM_uu * view(eigen_vck.vectors, 1:Nvck_uu, :)
    vectors_ijR_dd = BM_dd * view(eigen_vck.vectors, (Nvck_uu+1):(Nvck_uu+Nvck_dd), :)

    eigen_ijR = Eigen(copy(eigen_vck.values), [vectors_ijR_uu; vectors_ijR_dd])
    inv_sqrt_ndegen = map(n->1/sqrt(n), bse.ndegen)
    double_qR = map(x -> 2 * (q ⋅ x), bse.Rgrid)
    double_qorb_up = map(x -> 2 * (q ⋅ x), bse.orb_location_up)
    double_qorb_dn = map(x -> 2 * (q ⋅ x), bse.orb_location_dn)

    return BSEspinaware_wiTDA_uEigen_Sz0{typeof(eigen_ijR),qT}(
        q, eigen_ijR, inv_sqrt_ndegen, double_qR, double_qorb_up, double_qorb_dn, bse.ijR_uu, bse.ijR_dd)
end
struct BSEspinaware_wiTDA_uEigen_Sznp1{eigenT<:Eigen,qT<:Real} <: BSE_uEigen
    q::ReducedCoordinates{qT}
    eigen::eigenT
    inv_sqrt_ndegen::Vector{Float64}
    double_qR::Vector{qT}
    double_qorb_v::Vector{Float64} # up for n1, dn for p1
    double_qorb_c::Vector{Float64} # dn for n1, up for p1
    ijR::TripleList{Int,Int,Int} # ud for n1, du for p1
end
function (u::BSEspinaware_wiTDA_uEigen_Sznp1)(λ::Real, bandindex::AbstractVector{<:Integer}=eachindex(u.eigen.values))

    eqR = cispi.(-λ * u.double_qR) .* u.inv_sqrt_ndegen
    δ = 1 - λ
    eqorb₀ = [cispi(-λ * x - δ * y) for x in u.double_qorb_v, y in u.double_qorb_c]

    NijR = length(u.ijR)
    eigen_vectors = Matrix{eltype(u.eigen.vectors)}(undef, NijR, length(bandindex))
    for idx in Base.OneTo(NijR)
        (i, j, iR) = u.ijR[idx]
        eigen_vectors[idx, :] .= u.eigen.vectors[idx, bandindex] .* (eqorb₀[i, j] * eqR[iR])
    end

    return Eigen(u.eigen.values[bandindex], eigen_vectors)
end
function _BSE_uEigen(bse::BSEspinaware_wiTDA, q::ReducedCoordinates{qT}, eigen_vck::Eigen,
    bandkq_up, bandkq_dn, e_kR, ::Val{-1}) where {qT}
    return _BSE_uEigen(bse, q, eigen_vck, bandkq_dn, e_kR, Val(-1))
end
function _BSE_uEigen(bse::BSEspinaware_wiTDA, q::ReducedCoordinates{qT}, eigen_vck::Eigen,
    bandkq_dn, e_kR, ::Val{-1}) where {qT}

    BM = _vck2ijR(bse, bandkq_dn, e_kR, Val(-1))

    eigen_ijR = Eigen(deepcopy(eigen_vck.values), BM * eigen_vck.vectors)
    inv_sqrt_ndegen = map(n->1/sqrt(n), bse.ndegen)
    double_qR = map(x -> 2 * (q ⋅ x), bse.Rgrid)
    double_qorb_up = map(x -> 2 * (q ⋅ x), bse.orb_location_up)
    double_qorb_dn = map(x -> 2 * (q ⋅ x), bse.orb_location_dn)

    return BSEspinaware_wiTDA_uEigen_Sznp1{typeof(eigen_ijR),qT}(
        q, eigen_ijR, inv_sqrt_ndegen, double_qR, double_qorb_up, double_qorb_dn, bse.ijR_ud)
end
function _BSE_uEigen(bse::BSEspinaware_wiTDA, q::ReducedCoordinates{qT}, eigen_vck::Eigen,
    bandkq_up, bandkq_dn, e_kR, ::Val{1}) where {qT}
    return _BSE_uEigen(bse, q, eigen_vck, bandkq_up, e_kR, Val(1))
end
function _BSE_uEigen(bse::BSEspinaware_wiTDA, q::ReducedCoordinates{qT}, eigen_vck::Eigen,
    bandkq_up, e_kR, ::Val{1}) where {qT}

    BM = _vck2ijR(bse, bandkq_up, e_kR, Val(1))

    eigen_ijR = Eigen(deepcopy(eigen_vck.values), BM * eigen_vck.vectors)
    inv_sqrt_ndegen = map(n->1/sqrt(n), bse.ndegen)
    double_qR = map(x -> 2 * (q ⋅ x), bse.Rgrid)
    double_qorb_up = map(x -> 2 * (q ⋅ x), bse.orb_location_up)
    double_qorb_dn = map(x -> 2 * (q ⋅ x), bse.orb_location_dn)

    return BSEspinaware_wiTDA_uEigen_Sznp1{typeof(eigen_ijR),qT}(
        q, eigen_ijR, inv_sqrt_ndegen, double_qR, double_qorb_dn, double_qorb_up, bse.ijR_du)
end
