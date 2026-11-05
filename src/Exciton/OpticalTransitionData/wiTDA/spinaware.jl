function WannierInterActionBase.OpticalTransitionData(bse::BSEspinaware_wiTDA)

    ΔE = Vector{Float64}(undef, 0) # eV
    berryconnection = Matrix{ComplexF64}(undef, 3, 0) # Å

    Nvck_uu = length(bse.vck_uu)
    if Nvck_uu ≠ 0
        scissor_c_up = bse.scissor_c_up
        scissor_v_up = bse.scissor_v_up
        scissor_uu = scissor_c_up - scissor_v_up
        ΔE_uu = Vector{Float64}(undef, Nvck_uu)
        for i in Base.OneTo(Nvck_uu)
            (v, c, k) = bse.vck_uu[i]
            ΔE_uu[i] = bse.bandk_up[k].values[c] - bse.bandk_up[k].values[v] + scissor_uu
        end
        ΔE = [ΔE; ΔE_uu]
        berryconnection = [berryconnection bse.berryconnection_up]
    end
    Nvck_dd = length(bse.vck_dd)
    if Nvck_dd ≠ 0
        scissor_c_dn = bse.scissor_c_dn
        scissor_v_dn = bse.scissor_v_dn
        scissor_dd = scissor_c_dn - scissor_v_dn
        ΔE_dd = Vector{Float64}(undef, Nvck_dd)
        for i in Base.OneTo(Nvck_dd)
            (v, c, k) = bse.vck_dd[i]
            ΔE_dd[i] = bse.bandk_dn[k].values[c] - bse.bandk_dn[k].values[v] + scissor_dd
        end
        ΔE = [ΔE; ΔE_dd]
        berryconnection = [berryconnection bse.berryconnection_dn]
    end

    nk = length(bse.kgrid)
    Ω = abs(det(parent(bse.lattice)))
    NΩ = nk * Ω

    np = count(bse.period)
    if np == 3
        return OpticalTransitionData(ΔE, berryconnection; nspin=2, nspinor=1, NΩ)
    elseif np == 2
        period = bse.period
        if !period[1]
            Lz = bse.lattice[1, 1]
        elseif !period[2]
            Lz = bse.lattice[2, 2]
        else
            Lz = bse.lattice[3, 3]
        end
        return OpticalTransitionData(ΔE, berryconnection; nspin=2, nspinor=1, NΩ, Lz)
    end
end
function WannierInterActionBase.OpticalTransitionData(bse::BSEspinaware_wiTDA, eigen_vck::Eigen)

    ΔEA = Matrix{ComplexF64}(undef, 3, 0) # eV ⋅ Å

    Nvck_uu = length(bse.vck_uu)
    if Nvck_uu ≠ 0
        ΔEA_uu = Matrix{ComplexF64}(undef, 3, Nvck_uu)
        for i in Base.OneTo(Nvck_uu)
            (v, c, k) = bse.vck_uu[i]
            ΔE_vc = bse.bandk_up[k].values[c] - bse.bandk_up[k].values[v]
            ΔEA_uu[:, i] = bse.berryconnection_up[:, i] .* ΔE_vc
        end
        ΔEA = [ΔEA ΔEA_uu]
    end
    Nvck_dd = length(bse.vck_dd)
    if Nvck_dd ≠ 0
        ΔEA_dd = Matrix{ComplexF64}(undef, 3, Nvck_dd)
        for i in Base.OneTo(Nvck_dd)
            (v, c, k) = bse.vck_dd[i]
            ΔE_vc = bse.bandk_dn[k].values[c] - bse.bandk_dn[k].values[v]
            ΔEA_dd[:, i] = bse.berryconnection_dn[:, i] .* ΔE_vc
        end
        ΔEA = [ΔEA ΔEA_dd]
    end

    ΔE = copy(eigen_vck.values) # eV
    berryconnection = ΔEA * eigen_vck.vectors # eV ⋅ Å
    for i in eachindex(ΔE)
        berryconnection[:, i] ./= ΔE[i] # Å
    end

    nk = length(bse.kgrid)
    Ω = abs(det(parent(bse.lattice)))
    NΩ = nk * Ω

    np = count(bse.period)
    if np == 3
        return OpticalTransitionData(ΔE, berryconnection; nspin=2, nspinor=1, NΩ)
    elseif np == 2
        period = bse.period
        if !period[1]
            Lz = bse.lattice[1, 1]
        elseif !period[2]
            Lz = bse.lattice[2, 2]
        else
            Lz = bse.lattice[3, 3]
        end
        return OpticalTransitionData(ΔE, berryconnection; nspin=2, nspinor=1, NΩ, Lz)
    end
end
