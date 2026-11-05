function WannierInterActionBase.OpticalTransitionData(bse::BSEspinblind_wiTDA)

    Nvck = length(bse.vck)
    ΔE = Vector{Float64}(undef, Nvck) # eV
    for i in Base.OneTo(Nvck)
        (v, c, ik) = bse.vck[i]
        ΔE[i] = bse.bandk[ik].values[c] - bse.bandk[ik].values[v] + bse.scissor
    end

    nk = length(bse.kgrid)
    Ω = abs(det(parent(bse.lattice)))
    NΩ = nk * Ω

    if bse.channel == 1
        nspinor = 2 # for spinor or single spin.
    else
        nspinor = 1
    end

    np = count(bse.period)
    if np == 3
        return OpticalTransitionData(ΔE, bse.berryconnection; nspin=1, nspinor, NΩ)
    elseif np == 2
        period = bse.period
        if !period[1]
            Lz = bse.lattice[1, 1]
        elseif !period[2]
            Lz = bse.lattice[2, 2]
        else
            Lz = bse.lattice[3, 3]
        end
        return OpticalTransitionData(ΔE, bse.berryconnection; nspin=1, nspinor, NΩ, Lz)
    end
end
function WannierInterActionBase.OpticalTransitionData(bse::BSEspinblind_wiTDA, eigen_vck::Eigen)

    Nvck = length(bse.vck)

    ΔEA = Matrix{ComplexF64}(undef, 3, Nvck) # eV ⋅ Å
    for i in Base.OneTo(Nvck)
        (v, c, ik) = bse.vck[i]
        ΔE_vc = bse.bandk[ik].values[c] - bse.bandk[ik].values[v]
        ΔEA[:, i] = bse.berryconnection[:, i] .* ΔE_vc
    end

    ΔE = copy(eigen_vck.values) # eV
    berryconnection = ΔEA * eigen_vck.vectors # eV ⋅ Å
    for i in eachindex(ΔE)
        berryconnection[:, i] ./= ΔE[i] # Å
    end

    nk = length(bse.kgrid)
    Ω = abs(det(parent(bse.lattice)))
    NΩ = nk * Ω

    if bse.channel == 1
        nspinor = 2 # for spinor or single spin.
    else
        nspinor = 1
    end

    np = count(bse.period)
    if np == 3
        return OpticalTransitionData(ΔE, berryconnection; nspin=1, nspinor, NΩ)
    elseif np == 2
        period = bse.period
        if !period[1]
            Lz = bse.lattice[1, 1]
        elseif !period[2]
            Lz = bse.lattice[2, 2]
        else
            Lz = bse.lattice[3, 3]
        end
        return OpticalTransitionData(ΔE, berryconnection; nspin=1, nspinor, NΩ, Lz)
    end
end
