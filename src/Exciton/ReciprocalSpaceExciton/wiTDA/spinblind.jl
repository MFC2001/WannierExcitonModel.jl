struct ReciprocalSpaceBlochExciton_spinblind_wiTDA{eleHT<:HermitianReciprocalHoppings} <: ReciprocalSpaceExciton
    q::ReducedCoordinates{Float64}
    AijR::Matrix{ComplexF64} # [NijR, nband]
    ijR::TripleList{Int,Int,Int}
    norb::Int
    Rgrid::Vector{ReducedCoordinates{Int}}
    inv_ndegen::Vector{Float64}
    nv::Int
    nc::Int
    vc::Vector{Tuple{Int,Int}}
    eleH::eleHT
end
function (φk::ReciprocalSpaceBlochExciton_spinblind_wiTDA)(k::AbstractVector{<:Real},
    bandindex::Union{Integer,AbstractVector{<:Integer},Nothing}=nothing)
    if isnothing(bandindex)
        if size(φk.AijR, 2) == 1
            bandindex = 1
        else
            bandindex = axes(φk.AijR, 2)
        end
    end
    k = ReducedCoordinates(k)
    bandk = eigvecs(φk.eleH(k))
    bandkq = eigvecs(φk.eleH(k+φk.q))
    phase_normalize!(:positive_maxabs, bandk)
    phase_normalize!(:positive_maxabs, bandkq)
    e_kR = map(RR->cispi(2 * (k ⋅ RR)), φk.Rgrid)

    Nvc = length(φk.vc)
    NijR = length(φk.ijR)
    norb = φk.norb
    BM = Matrix{ComplexF64}(undef, NijR, Nvc)
    buffer_UU = Matrix{ComplexF64}(undef, norb, norb)
    for ivc in Base.OneTo(Nvc)
        (v, c) = φk.vc[ivc]
        ψv = view(bandk, :, v)
        ψc = view(bandkq, :, c)
        for j in Base.OneTo(norb), i in Base.OneTo(norb)
            buffer_UU[i, j] = ψv[i] * conj(ψc[j])
        end
        @simd for idx in Base.OneTo(NijR)
            (i, j, iR) = φk.ijR[idx]
            BM[idx, ivc] = buffer_UU[i, j] * e_kR[iR] * φk.inv_ndegen[iR]
        end
    end
    BM = permutedims(BM)
    AijR = view(φk.AijR, :, bandindex)
    Avc = BM * AijR
    Avc = reshape(Avc, φk.nv, φk.nc, :)
    if bandindex isa Integer
        return reshape(Avc, φk.nv, φk.nc)
    else
        return Avc
    end
end
function (bse::BSEspinblind_wiTDA)(::Val{:interpolation}, q::AbstractVector{<:Real}, Avck::AbstractVecOrMat{<:Number})
    size(Avck, 1) == length(bse.vck) || throw(DimensionMismatch(
        "vector length should be Nvck=$(length(bse.vck)), got $(size(Avck, 1))"))
    BM = bse(Val(:vck2ijR), q)
    AijR = BM * Avck
    if AijR isa AbstractVector
        AijR = reshape(AijR, :, 1)
    end
    inv_ndegen = 1 ./ bse.ndegen

    v = values(bse.vck)[1]
    c = values(bse.vck)[2]
    nv = length(v)
    nc = length(c)
    vc = Matrix{Tuple{Int,Int}}(undef, nv, nc)
    for ic in eachindex(c), iv in eachindex(v)
        vc[iv, ic] = (v[iv], c[ic])
    end
    vc = reshape(vc, :)

    return ReciprocalSpaceBlochExciton_spinblind_wiTDA(
        q, AijR, bse.ijR, numorb(bse), bse.Rgrid, inv_ndegen, nv, nc, vc, bse.eleH)
end
