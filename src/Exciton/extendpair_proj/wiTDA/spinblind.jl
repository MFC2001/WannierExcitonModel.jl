"""
```julia
julia> bse(Val(:extendpair_proj), q, wave, pairs...)
```

Project each eigenstate `wave[:, ib]` (coefficients in the `vck` basis) onto the
extended pair basis `|q, ijR>`, and return the summed weight

    P_{ijR}^{αq} = (1/√N) ∑_{vck} A_{vck}^{αq} U_{iv}^{k*} U_{jc}^{k+q} e^{i(k+q)·R},
    weight[ib]  = ∑_{(i,j,R) ∈ pairs} |P_{ijR}^{αq}|².

Each element of `pairs` is a tuple `(i, j, R)`, where `i` (`j`) is the hole
(electron) orbital and `R` is the electron-hole separation in reduced
coordinates.
"""
function (bse::BSEspinblind_wiTDA)(::Val{:extendpair_proj}, q::AbstractVector{<:Real},
    wave::AbstractVecOrMat{<:Number}, pairs::WanIntPairVecType...)
    kgrid = bse.kgrid
    bandk = bse.bandk
    vck = bse.vck
    Nvck = length(vck)
    size(wave, 1) == Nvck || throw(DimensionMismatch(
        "vector length should be Nvck=$Nvck, got $(size(wave, 1))"))

    q = ReducedCoordinates(q)
    bandkq = BAND(map(k -> k + q, kgrid), bse.eleH; vector=true)
    # keep the same gauge as the one used when building the BSE Hamiltonian
    phase_normalize!(:positive_maxabs, bandkq)

    np = length(pairs)
    nb = size(wave, 2)
    inv_sqrt_nk = 1 / sqrt(length(kgrid))
    P = Matrix{ComplexF64}(undef, nb, np)

    # coefficients multiplying A_{vck} in P_{ijR}; they do not depend on band
    coefs = Vector{ComplexF64}(undef, Nvck)
    for ip in Base.OneTo(np)
        (i, j, R) = pairs[ip]
        R = ReducedCoordinates{Int}(R)
        e_qR = cispi(2 * (q ⋅ R))
        e_kR = [cispi(2 * (k ⋅ R)) for k in kgrid]
        @inbounds for ivck in Base.OneTo(Nvck)
            (v, c, ik) = vck[ivck]
            coefs[ivck] = conj(bandk[ik].vectors[i, v]) * bandkq[ik].vectors[j, c] * (e_qR * e_kR[ik])
        end
        for ib in Base.OneTo(nb)
            Pib = zero(ComplexF64)
            @inbounds for ivck in Base.OneTo(Nvck)
                Pib += coefs[ivck] * wave[ivck, ib]
            end
            P[ib, ip] = inv_sqrt_nk * Pib
        end
    end

    # per-state projection weight, ∑_{ijR∈F} |P_{ijR}|²
    return vec(sum(abs2, P; dims=2))
end
