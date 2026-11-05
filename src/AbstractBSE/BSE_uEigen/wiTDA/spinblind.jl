"""
```julia
julia> H = bse(q)
julia> u = bse(:periodic, q, eigen!(H))
```
"""
function (bse::BSEspinblind_wiTDA)(::Val{:periodic}, q::AbstractVector{<:Real}, eigen_vck::Eigen)
    size(eigen_vck.vectors, 1) == length(bse.vck) || throw(DimensionMismatch(
        "vector length should be Nvck=$(length(bse.vck)), got $(size(eigen_vck.vectors, 1))"))
    # bandkq
    q = ReducedCoordinates(q)
    bandkq = BAND(map(k -> k + q, bse.kgrid), bse.eleH; vector=true)
    phase_normalize!(:positive_maxabs, bandkq)
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
    return _BSE_uEigen(bse, q, eigen_vck, bandkq, e_kR)
end
struct BSEspinblind_wiTDA_uEigen{eigenT<:Eigen,qT<:Real} <: BSE_uEigen
    q::ReducedCoordinates{qT}
    eigen::eigenT
    inv_sqrt_ndegen::Vector{Float64}
    double_qR::Vector{qT}
    double_qorb::Vector{Float64}
    ijR::TripleList{Int,Int,Int}
end
function (u::BSEspinblind_wiTDA_uEigen)(λ::Real, bandindex::AbstractVector{<:Integer}=eachindex(u.eigen.values))

    eqR = cispi.((-λ) .* u.double_qR) .* u.inv_sqrt_ndegen
    δ = 1 - λ
    eqorb₀ = [cispi(-λ * x - δ * y) for x in u.double_qorb, y in u.double_qorb]

    NijR = length(u.ijR)
    eigen_vectors = Matrix{eltype(u.eigen.vectors)}(undef, NijR, length(bandindex))
    Threads.@threads for idx in Base.OneTo(NijR)
        (i, j, iR) = u.ijR[idx]
        eigen_vectors[idx, :] .= u.eigen.vectors[idx, bandindex] .* (eqorb₀[i, j] * eqR[iR])
    end

    return Eigen(u.eigen.values[bandindex], eigen_vectors)
end
function _BSE_uEigen(bse::BSEspinblind_wiTDA, q::ReducedCoordinates{qT}, eigen_vck::Eigen, bandkq, e_kR) where {qT}
    BM = _vck2ijR(bse, bandkq, e_kR)
    eigen_ijR = Eigen(copy(eigen_vck.values), BM * eigen_vck.vectors)
    inv_sqrt_ndegen = map(n->1/sqrt(n), bse.ndegen)
    double_qR = map(x -> 2 * (q ⋅ x), bse.Rgrid)
    double_qorb = map(x -> 2 * (q ⋅ x), bse.orb_location)
    return BSEspinblind_wiTDA_uEigen{typeof(eigen_ijR),qT}(
        q, eigen_ijR, inv_sqrt_ndegen, double_qR, double_qorb, bse.ijR)
end
