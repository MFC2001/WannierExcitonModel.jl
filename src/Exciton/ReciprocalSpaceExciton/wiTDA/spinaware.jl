"""
    (bse::BSEspinaware_wiTDA)(:interpolation, q::AbstractVector{<:Real}, Avck::AbstractVecOrMat{<:Number}, Sz::Integer)

Return a callable object that computes the expansion coefficients, in the band basis, of the exciton state(s) at an arbitrary momentum point in reciprocal space.

`q` is the exciton momentum corresponding to `Avck`. 
`Avck` is the exciton wavefunction; multiple exciton states may be passed at once, in which case its second dimension is the number of exciton states.
`Sz` is the z-component of the exciton's total spin.

# Example
```julia
julia> φ = bse(:interpolation, q, Avck, 0);
julia> (Avc_uu, Avc_dd) = φ([0, 0, 0]);
julia> size(Avc_uu) == (nv_up, nc_up, size(Avck, 2))

true
julia> φ = bse(:interpolation, q, Avck, 1);
julia> size(φ([0, 0, 0])) == (nv_dn, nc_up, size(Avck, 2))

true
```
"""
function (bse::BSEspinaware_wiTDA)(::Val{:interpolation}, q::AbstractVector{<:Real}, Avck::AbstractVecOrMat{<:Number}, Sz::Integer)
    q = ReducedCoordinates(q)
    return _ReciprocalSpaceBlochExciton_spinaware_wiTDA(bse, q, Avck, Val(Sz))
end
struct ReciprocalSpaceBlochExciton_Sz0_wiTDA{eleHT_up,eleHT_dn} <: ReciprocalSpaceExciton
    q::ReducedCoordinates{Float64}
    AijR_uu::Matrix{ComplexF64} # [NijR_uu, nband]
    AijR_dd::Matrix{ComplexF64} # [NijR_dd, nband]
    ijR_uu::TripleList{Int,Int,Int}
    ijR_dd::TripleList{Int,Int,Int}
    norb_up::Int
    norb_dn::Int
    Rgrid::Vector{ReducedCoordinates{Int}}
    inv_ndegen::Vector{Float64}
    nv_up::Int
    nc_up::Int
    vc_uu::Vector{Tuple{Int,Int}}
    nv_dn::Int
    nc_dn::Int
    vc_dd::Vector{Tuple{Int,Int}}
    eleH_up::eleHT_up
    eleH_dn::eleHT_dn
end
function (φk::ReciprocalSpaceBlochExciton_Sz0_wiTDA)(k::AbstractVector{<:Real},
    bandindex::Union{Integer,AbstractVector{<:Integer},Nothing}=nothing)
    if isnothing(bandindex)
        if size(φk.AijR_uu, 2) == 1
            bandindex = 1
        else
            bandindex = axes(φk.AijR_uu, 2)
        end
    end
    k = ReducedCoordinates(k)
    e_kR = map(RR->cispi(2 * (k ⋅ RR)), φk.Rgrid)
    if !isempty(φk.AijR_uu)
        bandk_up = eigvecs(φk.eleH_up(k))
        bandkq_up = eigvecs(φk.eleH_up(k+φk.q))
        phase_normalize!(:positive_maxabs, bandk_up)
        phase_normalize!(:positive_maxabs, bandkq_up)
        Nvc_uu = length(φk.vc_uu)
        NijR_uu = length(φk.ijR_uu)
        norb_up = φk.norb_up
        BM = Matrix{ComplexF64}(undef, NijR_uu, Nvc_uu)
        buffer_UU = Matrix{ComplexF64}(undef, norb_up, norb_up)
        for ivc in Base.OneTo(Nvc_uu)
            (v, c) = φk.vc_uu[ivc]
            ψv = view(bandk_up, :, v)
            ψc = view(bandkq_up, :, c)
            for j in Base.OneTo(norb_up), i in Base.OneTo(norb_up)
                buffer_UU[i, j] = ψv[i] * conj(ψc[j])
            end
            @simd for idx in Base.OneTo(NijR_uu)
                (i, j, iR) = φk.ijR_uu[idx]
                BM[idx, ivc] = buffer_UU[i, j] * e_kR[iR] * φk.inv_ndegen[iR]
            end
        end
        BM = permutedims(BM)
        AijR_uu = view(φk.AijR_uu, :, bandindex)
        Avc_uu = BM * AijR_uu
        Avc_uu = reshape(Avc_uu, φk.nv_up, φk.nc_up, :)
    else
        # it should be empty
        Avc_uu = Array{ComplexF64}(undef, φk.nv_up, φk.nc_up, length(bandindex))
    end
    if bandindex isa Integer
        Avc_uu = reshape(Avc_uu, φk.nv_up, φk.nc_up)
    end
    if !isempty(φk.AijR_dd)
        bandk_dn = eigvecs(φk.eleH_dn(k))
        bandkq_dn = eigvecs(φk.eleH_dn(k+φk.q))
        phase_normalize!(:positive_maxabs, bandk_dn)
        phase_normalize!(:positive_maxabs, bandkq_dn)
        Nvc_dd = length(φk.vc_dd)
        NijR_dd = length(φk.ijR_dd)
        norb_dn = φk.norb_dn
        BM = Matrix{ComplexF64}(undef, NijR_dd, Nvc_dd)
        buffer_UU = Matrix{ComplexF64}(undef, norb_dn, norb_dn)
        for ivc in Base.OneTo(Nvc_dd)
            (v, c) = φk.vc_dd[ivc]
            ψv = view(bandk_dn, :, v)
            ψc = view(bandkq_dn, :, c)
            for j in Base.OneTo(norb_dn), i in Base.OneTo(norb_dn)
                buffer_UU[i, j] = ψv[i] * conj(ψc[j])
            end
            @simd for idx in Base.OneTo(NijR_dd)
                (i, j, iR) = φk.ijR_dd[idx]
                BM[idx, ivc] = buffer_UU[i, j] * e_kR[iR] * φk.inv_ndegen[iR]
            end
        end
        BM = permutedims(BM)
        AijR_dd = view(φk.AijR_dd, :, bandindex)
        Avc_dd = BM * AijR_dd
        Avc_dd = reshape(Avc_dd, φk.nv_dn, φk.nc_dn, :)
    else
        # it should be empty
        Avc_dd = Array{ComplexF64}(undef, φk.nv_dn, φk.nc_dn, length(bandindex))
    end
    if bandindex isa Integer
        Avc_dd = reshape(Avc_dd, φk.nv_dn, φk.nc_dn)
    end
    return Avc_uu, Avc_dd
end
function _ReciprocalSpaceBlochExciton_spinaware_wiTDA(bse::BSEspinaware_wiTDA, q::ReducedCoordinates{<:Real}, Avck::AbstractVecOrMat{<:Number}, ::Val{0})
    Nvck_uu = length(bse.vck_uu)
    Nvck_dd = length(bse.vck_dd)
    size(Avck, 1) == Nvck_uu+Nvck_dd || throw(DimensionMismatch(
        "length of exciton wavefunction should be Nvck=$(Nvck_uu+Nvck_dd), got $(size(Avck, 1))"))

    BM_uu, BM_dd = bse(Val(:vck2ijR), q, 0)
    if iszero(Nvck_uu)
        AijR_uu = Matrix{ComplexF64}(undef, 0, size(Avck, 2))
        nv_up = 0
        nc_up = 0
        vc_uu = Vector{Tuple{Int,Int}}(undef, 0)
    else
        AijR_uu = BM_uu * view(Avck, 1:Nvck_uu, :)
        if AijR_uu isa AbstractVector
            AijR_uu = reshape(AijR_uu, :, 1)
        end
        v_up = values(bse.vck_uu)[1]
        c_up = values(bse.vck_uu)[2]
        nv_up = length(v_up)
        nc_up = length(c_up)
        vc_uu = Matrix{Tuple{Int,Int}}(undef, nv_up, nc_up)
        for ic in eachindex(c_up), iv in eachindex(v_up)
            vc_uu[iv, ic] = (v_up[iv], c_up[ic])
        end
        vc_uu = reshape(vc_uu, :)
    end
    if iszero(Nvck_dd)
        AijR_dd = Matrix{ComplexF64}(undef, 0, size(Avck, 2))
        nv_dn = 0
        nc_dn = 0
        vc_dd = Vector{Tuple{Int,Int}}(undef, 0)
    else
        AijR_dd = BM_dd * view(Avck, (Nvck_uu+1):(Nvck_uu+Nvck_dd), :)
        if AijR_dd isa AbstractVector
            AijR_dd = reshape(AijR_dd, :, 1)
        end
        v_dn = values(bse.vck_dd)[1]
        c_dn = values(bse.vck_dd)[2]
        nv_dn = length(v_dn)
        nc_dn = length(c_dn)
        vc_dd = Matrix{Tuple{Int,Int}}(undef, nv_dn, nc_dn)
        for ic in eachindex(c_dn), iv in eachindex(v_dn)
            vc_dd[iv, ic] = (v_dn[iv], c_dn[ic])
        end
        vc_dd = reshape(vc_dd, :)
    end
    inv_ndegen = 1 ./ bse.ndegen

    return ReciprocalSpaceBlochExciton_Sz0_wiTDA{typeof(bse.eleH_up),typeof(bse.eleH_dn)}(
        q, AijR_uu, AijR_dd, bse.ijR_uu, bse.ijR_dd, numorb(bse.eleH_up), numorb(bse.eleH_dn),
        bse.Rgrid, inv_ndegen, nv_up, nc_up, vc_uu, nv_dn, nc_dn, vc_dd, bse.eleH_up, bse.eleH_dn)
end
struct ReciprocalSpaceBlochExciton_Sznp1_wiTDA{eleHT_v,eleHT_c} <: ReciprocalSpaceExciton
    q::ReducedCoordinates{Float64}
    AijR::Matrix{ComplexF64} # [NijR, nband]
    ijR::TripleList{Int,Int,Int}
    norb_v::Int
    norb_c::Int
    Rgrid::Vector{ReducedCoordinates{Int}}
    inv_ndegen::Vector{Float64}
    nv::Int
    nc::Int
    vc::Vector{Tuple{Int,Int}}
    eleH_v::eleHT_v
    eleH_c::eleHT_c
end
function (φk::ReciprocalSpaceBlochExciton_Sznp1_wiTDA)(k::AbstractVector{<:Real},
    bandindex::Union{Integer,AbstractVector{<:Integer},Nothing}=nothing)
    if isnothing(bandindex)
        if size(φk.AijR, 2) == 1
            bandindex = 1
        else
            bandindex = axes(φk.AijR, 2)
        end
    end
    k = ReducedCoordinates(k)
    e_kR = map(RR->cispi(2 * (k ⋅ RR)), φk.Rgrid)
    bandk = eigvecs(φk.eleH_v(k))
    bandkq = eigvecs(φk.eleH_c(k+φk.q))
    phase_normalize!(:positive_maxabs, bandk)
    phase_normalize!(:positive_maxabs, bandkq)
    Nvc = length(φk.vc)
    NijR = length(φk.ijR)
    norb_v = φk.norb_v
    norb_c = φk.norb_c
    BM = Matrix{ComplexF64}(undef, NijR, Nvc)
    buffer_UU = Matrix{ComplexF64}(undef, norb_v, norb_c)
    for ivc in Base.OneTo(Nvc)
        (v, c) = φk.vc[ivc]
        ψv = view(bandk, :, v)
        ψc = view(bandkq, :, c)
        for j in Base.OneTo(norb_c), i in Base.OneTo(norb_v)
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
function _ReciprocalSpaceBlochExciton_spinaware_wiTDA(bse::BSEspinaware_wiTDA, q::ReducedCoordinates{<:Real}, Avck::AbstractVecOrMat{<:Number}, ::Val{-1})
    size(Avck, 1) == length(bse.vck_ud) || throw(DimensionMismatch(
        "length of exciton wavefunction should be Nvck=$(length(bse.vck_ud)), got $(size(Avck, 1))"))

    BM = bse(Val(:vck2ijR), q, -1)
    AijR = BM * Avck
    if AijR isa AbstractVector
        AijR = reshape(AijR, :, 1)
    end
    v = values(bse.vck_ud)[1]
    c = values(bse.vck_ud)[2]
    nv = length(v)
    nc = length(c)
    vc = Matrix{Tuple{Int,Int}}(undef, nv, nc)
    for ic in eachindex(c), iv in eachindex(v)
        vc[iv, ic] = (v[iv], c[ic])
    end
    vc = reshape(vc, :)
    inv_ndegen = 1 ./ bse.ndegen

    return ReciprocalSpaceBlochExciton_Sznp1_wiTDA{typeof(bse.eleH_up),typeof(bse.eleH_dn)}(
        q, AijR, bse.ijR_ud, numorb(bse.eleH_up), numorb(bse.eleH_dn),
        bse.Rgrid, inv_ndegen, nv, nc, vc, bse.eleH_up, bse.eleH_dn)
end
function _ReciprocalSpaceBlochExciton_spinaware_wiTDA(bse::BSEspinaware_wiTDA, q::ReducedCoordinates{<:Real}, Avck::AbstractVecOrMat{<:Number}, ::Val{1})
    size(Avck, 1) == length(bse.vck_du) || throw(DimensionMismatch(
        "length of eaciton wavefunction should be Nvck=$(length(bse.vck_du)), got $(size(Avck, 1))"))

    BM = bse(Val(:vck2ijR), q, 1)
    AijR = BM * Avck
    if AijR isa AbstractVector
        AijR = reshape(AijR, :, 1)
    end
    v = values(bse.vck_du)[1]
    c = values(bse.vck_du)[2]
    nv = length(v)
    nc = length(c)
    vc = Matrix{Tuple{Int,Int}}(undef, nv, nc)
    for ic in eachindex(c), iv in eachindex(v)
        vc[iv, ic] = (v[iv], c[ic])
    end
    vc = reshape(vc, :)
    inv_ndegen = 1 ./ bse.ndegen

    return ReciprocalSpaceBlochExciton_Sznp1_wiTDA{typeof(bse.eleH_dn),typeof(bse.eleH_up)}(
        q, AijR, bse.ijR_du, numorb(bse.eleH_dn), numorb(bse.eleH_up),
        bse.Rgrid, inv_ndegen, nv, nc, vc, bse.eleH_dn, bse.eleH_up)
end
