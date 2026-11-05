# function _BSE_pH_kgrid(H, kgrid, lattice, orb_location)

#     norb = length(orb_location)
#     dorb = Matrix{eltype(orb_location)}(undef, norb, norb)
#     for j in 1:norb, i in 1:norb
#         dorb[i, j] = orb_location[i] - orb_location[j]
#     end

#     pH = Array{ComplexF64}(undef, norb, norb, 3, length(kgrid))
#     Threads.@threads for ik in eachindex(kgrid)
#         k = kgrid[ik]
#         pHk = view(pH,:,:,:,ik)
#         H(Val(:partial), pHk, lattice, k, orb_location)
#         # later use ψk, so remove the phase factor here.
#         for j in 2:norb, i in 1:(j-1)
#             phase = cispi(2 * (k ⋅ dorb[i, j]))
#             pHk[i, j, :] .*= phase
#             pHk[j, i, :] .*= conj(phase)
#         end
#     end

#     return pH
# end
struct _BSE_Γdata
    berryconnection::Matrix{ComplexF64} # [1:3, ivck], Å
    CoulombScaledivnk::Float64
    rlattice::ReciprocalLattice{Float64}
end
function _BSE_berryconnection_vck(H, kgrid, bandk, lattice, orblocat, vckmap)
    if Nvck == 0
        return Matrix{ComplexF64}(undef, 3, 0)
    end

    pH = Array{ComplexF64}(undef, norb, norb, 3, length(kgrid))
    Threads.@threads for ik in eachindex(kgrid)
        k = kgrid[ik]
        pHk = view(pH,:,:,:,ik)
        H(Val(:partial), pHk, k, lattice)
    end
    ipH = map(x->complex(-imag(x), real(x)), pH)
    orblocat = map(τ->lattice * τ, orblocat)

    Nvck = length(vckmap)
    berrconnection = Array{ComplexF64}(undef, 3, Nvck) # Å
    Threads.@threads for i in Base.OneTo(Nvck)
        (v, c, k) = vckmap[i]
        ψv = view(bandk[k].vectors, :, v)
        ψc = view(bandk[k].vectors, :, c)
        ΔE = bandk[k].values[c] - bandk[k].values[v]
        for α in 1:3
            ipHk = view(ipH,:,:,α,k)
            berrconnection[α, i] = ψv ⋅ (ipHk * ψc) / ΔE
        end
        berrconnection[α, i] = sum(i -> conj(ψv[i]) * ψc[i] * orblocat[i], 1:norb)
    end

    return berrconnection


    # rlattice = reciprocal(TB.lattice)
    # (blist, weights) = kgrid_shell_blist(kgrid, rlattice; nsupercell, shell_tol, orth_tol)

    # adjust the length of blist, adjust weights simultaneously.
    # blist_norm = map(b->norm(rlattice * b), blist)
    # blist = map(eachindex(blist)) do ib
    #     blist[ib] / blist_norm[ib] * normb
    # end
    # weights = weights .* blist_norm .^ 2 ./ normb^2








    # Ω = abs(det(parent(TB.lattice)))
    # CoulombScaledivnk = CoulombScale * 4π / Ω / length(kgrid)

    # rlattice = reciprocal(TB.lattice)


    # if Nvck == 0
    #     upHu = Matrix{ComplexF64}(undef, 3, 0)
    #     upu = Matrix{ComplexF64}(undef, 3, 0)
    #     return _BSE_Γdata(upHu, upu, CoulombScaledivnk, rlattice)
    # end

    # pH = _BSE_pH_kgrid(TB.H, kgrid, TB.lattice, TB.orb_location)

    # upHu = Matrix{ComplexF64}(undef, 3, Nvck)
    # Threads.@threads for i in Base.OneTo(Nvck)
    #     (v, c, k) = vckmap[i]
    #     ψv = view(bandk[k].vectors, :, v)
    #     ψc = view(bandk[k].vectors, :, c)
    #     for α in 1:3
    #         pHk = view(pH,:,:,α,k)
    #         upHu[α, i] = ψv ⋅ (pHk * ψc) # eV ⋅ Å
    #     end
    # end

    # if count(TB.period) == 3
    #     upu = Matrix{ComplexF64}(undef, 3, Nvck) # Å
    #     Threads.@threads for i in Base.OneTo(Nvck)
    #         (v, c, k) = vckmap[i]
    #         ΔE = bandk[k].values[c] - bandk[k].values[v] # eV
    #         upu[:, i] .= upHu[:, i] ./ ΔE
    #     end
    # else
    #     upu = Matrix{ComplexF64}(undef, 3, 0)
    # end

    # return _BSE_Γdata(upHu, upu, CoulombScaledivnk, rlattice)
end
