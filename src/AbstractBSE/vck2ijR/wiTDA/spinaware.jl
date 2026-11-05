function (bse::BSEspinaware_wiTDA)(::Val{:vck2ijR}, q::AbstractVector{<:Real}, Sz::Integer)
    if Sz == 0
        Nvck_Sz0 = length(bse.vck_uu) + length(bse.vck_dd)
        if iszero(Nvck_Sz0)
            return Matrix{ComplexF64}(undef, 0, 0)
        end
    elseif Sz == -1
        Nvck_Szn1 = length(bse.vck_ud)
        if iszero(Nvck_Szn1)
            return Matrix{ComplexF64}(undef, 0, 0)
        end
    elseif Sz == 1
        Nvck_Szp1 = length(bse.vck_du)
        if iszero(Nvck_Szp1)
            return Matrix{ComplexF64}(undef, 0, 0)
        end
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
        return _vck2ijR(bse, bandkq_up, bandkq_dn, e_kR, Val(0))
    elseif Sz == -1
        bandkq_up = nothing
        bandkq_dn = BAND(kq, bse.eleH_dn; vector=true)
        phase_normalize!(:positive_maxabs, bandkq_dn)
        return _vck2ijR(bse, bandkq_dn, e_kR, Val(-1))
    else # Sz == 1
        bandkq_dn = nothing
        bandkq_up = BAND(kq, bse.eleH_up; vector=true)
        phase_normalize!(:positive_maxabs, bandkq_up)
        return _vck2ijR(bse, bandkq_up, e_kR, Val(1))
    end
end
function _vck2ijR(bse::BSEspinaware_wiTDA,
    bandkq_up::Union{Vector{<:Eigen},Vector{<:Matrix{<:Number}},Nothing},
    bandkq_dn::Union{Vector{<:Eigen},Vector{<:Matrix{<:Number}},Nothing},
    e_kR::Matrix{<:Number}, ::Val{0})

    if bandkq_up isa Vector{<:Eigen}
        bandkq_up = map(x->x.vectors, bandkq_up)
    end
    if bandkq_dn isa Vector{<:Eigen}
        bandkq_dn = map(x->x.vectors, bandkq_dn)
    end

    nk = length(bse.kgrid)
    inversesqrtnk = 1 / sqrt(nk)

    # uu
    Nvck_uu = length(bse.vck_uu)
    NijR_uu = length(bse.ijR_uu)
    if iszero(Nvck_uu)
        BM_uu = Matrix{ComplexF64}(undef, NijR_uu, 0)
    else
        norb_up = length(bse.orb_location_up)
        bandk_up = bse.bandk_up
        BM_uu = Matrix{ComplexF64}(undef, NijR_uu, Nvck_uu)
        ivck_chunk = _split_tasks(Nvck_uu, Threads.nthreads())
        Threads.@threads for ivcks in ivck_chunk
            buffer_UU = Matrix{ComplexF64}(undef, norb_up, norb_up)
            @inbounds for ivck in ivcks
                (v, c, ik) = bse.vck_uu[ivck]
                ψv = view(bandk_up[ik].vectors, :, v)
                ψc = view(bandkq_up[ik], :, c)
                for j in Base.OneTo(norb_up), i in Base.OneTo(norb_up)
                    buffer_UU[i, j] = inversesqrtnk * conj(ψv[i]) * ψc[j]
                end
                @simd for idx in Base.OneTo(NijR_uu)
                    (i, j, iR) = bse.ijR_uu[idx]
                    BM_uu[idx, ivck] = buffer_UU[i, j] * e_kR[iR, ik]
                end
            end
        end
    end

    # dd
    Nvck_dd = length(bse.vck_dd)
    NijR_dd = length(bse.ijR_dd)
    if iszero(Nvck_dd)
        BM_dd = Matrix{ComplexF64}(undef, NijR_dd, 0)
    else
        norb_dn = length(bse.orb_location_dn)
        bandk_dn = bse.bandk_dn
        BM_dd = Matrix{ComplexF64}(undef, NijR_dd, Nvck_dd)
        ivck_chunk = _split_tasks(Nvck_dd, Threads.nthreads())
        Threads.@threads for ivcks in ivck_chunk
            buffer_UU = Matrix{ComplexF64}(undef, norb_dn, norb_dn)
            @inbounds for ivck in ivcks
                (v, c, ik) = bse.vck_dd[ivck]
                ψv = view(bandk_dn[ik].vectors, :, v)
                ψc = view(bandkq_dn[ik], :, c)
                for j in Base.OneTo(norb_dn), i in Base.OneTo(norb_dn)
                    buffer_UU[i, j] = inversesqrtnk * conj(ψv[i]) * ψc[j]
                end
                @simd for idx in Base.OneTo(NijR_dd)
                    (i, j, iR) = bse.ijR_dd[idx]
                    BM_dd[idx, ivck] = buffer_UU[i, j] * e_kR[iR, ik]
                end
            end
        end
    end

    return BM_uu, BM_dd
end
function _vck2ijR(bse::BSEspinaware_wiTDA,
    bandkq_up::Union{Vector{<:Eigen},Vector{<:Matrix{<:Number}},Nothing},
    bandkq_dn::Union{Vector{<:Eigen},Vector{<:Matrix{<:Number}}},
    e_kR::Matrix{<:Number}, ::Val{-1})
    return _vck2ijR(bse, bandkq_dn, e_kR, Val(-1))
end
function _vck2ijR(bse::BSEspinaware_wiTDA,
    bandkq_dn::Union{Vector{<:Eigen},Vector{<:Matrix{<:Number}}}, e_kR::Matrix{<:Number}, ::Val{-1})

    if bandkq_dn isa Vector{<:Eigen}
        bandkq_dn = map(x->x.vectors, bandkq_dn)
    end

    nk = length(bse.kgrid)
    inversesqrtnk = 1 / sqrt(nk)
    bandk_up = bse.bandk_up

    # ud
    norb_up = length(bse.orb_location_up)
    norb_dn = length(bse.orb_location_dn)
    Nvck_ud = length(bse.vck_ud)
    NijR_ud = length(bse.ijR_ud)
    BM_ud = Matrix{ComplexF64}(undef, NijR_ud, Nvck_ud)
    ivck_chunk = _split_tasks(Nvck_ud, Threads.nthreads())
    Threads.@threads for ivcks in ivck_chunk
        buffer_UU = Matrix{ComplexF64}(undef, norb_up, norb_dn)
        @inbounds for ivck in ivcks
            (v, c, ik) = bse.vck_ud[ivck]
            ψv = view(bandk_up[ik].vectors, :, v)
            ψc = view(bandkq_dn[ik], :, c)
            for j in Base.OneTo(norb_dn), i in Base.OneTo(norb_up)
                buffer_UU[i, j] = inversesqrtnk * conj(ψv[i]) * ψc[j]
            end
            @simd for idx in Base.OneTo(NijR_ud)
                (i, j, iR) = bse.ijR_ud[idx]
                BM_ud[idx, ivck] = buffer_UU[i, j] * e_kR[iR, ik]
            end
        end
    end

    return BM_ud
end
function _vck2ijR(bse::BSEspinaware_wiTDA,
    bandkq_up::Union{Vector{<:Eigen},Vector{<:Matrix{<:Number}}},
    bandkq_dn::Union{Vector{<:Eigen},Vector{<:Matrix{<:Number}},Nothing},
    e_kR::Matrix{<:Number}, ::Val{1})
    return _vck2ijR(bse, bandkq_up, e_kR, Val(1))
end
function _vck2ijR(bse::BSEspinaware_wiTDA,
    bandkq_up::Union{Vector{<:Eigen},Vector{<:Matrix{<:Number}}}, e_kR::Matrix{<:Number}, ::Val{1})

    if bandkq_up isa Vector{<:Eigen}
        bandkq_up = map(x->x.vectors, bandkq_up)
    end

    nk = length(bse.kgrid)
    inversesqrtnk = 1 / sqrt(nk)
    bandk_dn = bse.bandk_dn

    # du
    norb_up = length(bse.orb_location_up)
    norb_dn = length(bse.orb_location_dn)
    Nvck_du = length(bse.vck_du)
    NijR_du = length(bse.ijR_du)
    BM_du = Matrix{ComplexF64}(undef, NijR_du, Nvck_du)
    ivck_chunk = _split_tasks(Nvck_du, Threads.nthreads())
    Threads.@threads for ivcks in ivck_chunk
        buffer_UU = Matrix{ComplexF64}(undef, norb_dn, norb_up)
        @inbounds for ivck in ivcks
            (v, c, ik) = bse.vck_ud[ivck]
            ψv = view(bandk_dn[ik].vectors, :, v)
            ψc = view(bandkq_up[ik], :, c)
            for j in Base.OneTo(norb_up), i in Base.OneTo(norb_dn)
                buffer_UU[i, j] = inversesqrtnk * conj(ψv[i]) * ψc[j]
            end
            @simd for idx in Base.OneTo(NijR_du)
                (i, j, iR) = bse.ijR_du[idx]
                BM_du[idx, ivck] = buffer_UU[i, j] * e_kR[iR, ik]
            end
        end
    end

    return BM_du
end
