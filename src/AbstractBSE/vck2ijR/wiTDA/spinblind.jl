function (bse::BSEspinblind_wiTDA)(::Val{:vck2ijR}, q::AbstractVector{<:Real})
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
    return _vck2ijR(bse, bandkq, e_kR)
end
function _vck2ijR(bse::BSEspinblind_wiTDA, bandkq::Vector{<:Eigen}, e_kR::Matrix{<:Number})
    return _vck2ijR(bse, map(x->x.vectors, bandkq), e_kR)
end
function _vck2ijR(bse::BSEspinblind_wiTDA, bandkq::Vector{<:Matrix{<:Number}}, e_kR::Matrix{<:Number})

    norb = numorb(bse)
    Nvck = length(bse.vck)
    NijR = length(bse.ijR)
    nk = length(bse.kgrid)
    inversesqrtnk = 1 / sqrt(nk)
    bandk = bse.bandk

    BM = Matrix{ComplexF64}(undef, NijR, Nvck)
    ivck_chunk = _split_tasks(Nvck, Threads.nthreads())
    Threads.@threads for ivcks in ivck_chunk
        buffer_UU = Matrix{ComplexF64}(undef, norb, norb)
        @inbounds for ivck in ivcks
            (v, c, ik) = bse.vck[ivck]
            ψv = view(bandk[ik].vectors, :, v)
            ψc = view(bandkq[ik], :, c)
            for j in Base.OneTo(norb), i in Base.OneTo(norb)
                buffer_UU[i, j] = inversesqrtnk * conj(ψv[i]) * ψc[j]
            end
            @simd for idx in Base.OneTo(NijR)
                (i, j, iR) = bse.ijR[idx]
                BM[idx, ivck] = buffer_UU[i, j] * e_kR[iR, ik]
            end
        end
    end

    return BM
end
