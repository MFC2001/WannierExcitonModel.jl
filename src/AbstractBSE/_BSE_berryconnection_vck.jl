function _BSE_berryconnection_vck(H, kgrid, bandk, lattice, orblocat, vckmap)
    Nvck = length(vckmap)
    if Nvck == 0
        return Matrix{ComplexF64}(undef, 3, 0)
    end

    norb = length(orblocat)
    pH = Array{ComplexF64}(undef, norb, norb, 3, length(kgrid))
    Threads.@threads for ik in eachindex(kgrid)
        k = kgrid[ik]
        pHk = view(pH, :, :, :, ik)
        H(Val(:partial), pHk, k, lattice)
    end
    ipH = map(x->complex(-imag(x), real(x)), pH)
    orblocat = map(τ->lattice * τ, orblocat)

    berrconnection = Array{ComplexF64}(undef, 3, Nvck) # Å
    Threads.@threads for ivck in Base.OneTo(Nvck)
        (v, c, k) = vckmap[ivck]
        ψv = view(bandk[k].vectors, :, v)
        ψc = view(bandk[k].vectors, :, c)
        ΔE = bandk[k].values[c] - bandk[k].values[v]
        for α in 1:3
            ipHαk = view(ipH, :, :, α, k)
            berrconnection[α, ivck] = ψv ⋅ (ipHαk * ψc) / ΔE
        end
        UUτ = zeros(ComplexF64, 3)
        for i in 1:norb
            UUτ += conj(ψv[i]) * ψc[i] * orblocat[i]
        end
        berrconnection[:, ivck] = berrconnection[:, ivck] + UUτ
    end

    return berrconnection
end
