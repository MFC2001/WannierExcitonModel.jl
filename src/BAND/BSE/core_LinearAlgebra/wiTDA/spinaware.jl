function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEspinaware_wiTDA, ::Val{false}, ::Val{:Bloch})

    Nq = length(qpoints)
    (Nvck_Sz0, Nvck_Szn1, Nvck_Szp1) = bse(Val(:dimension))
    has_Sz0 = !iszero(Nvck_Sz0)
    has_Szn1 = !iszero(Nvck_Szn1)
    has_Szp1 = !iszero(Nvck_Szp1)
    Sz = Int[]
    has_Sz0 && push!(Sz, 0)
    has_Szn1 && push!(Sz, -1)
    has_Szp1 && push!(Sz, 1)
    HSz0 = has_Sz0 ? Matrix{ComplexF64}(undef, Nvck_Sz0, Nvck_Sz0) : nothing
    HSzn1 = has_Szn1 ? Matrix{ComplexF64}(undef, Nvck_Szn1, Nvck_Szn1) : nothing
    HSzp1 = has_Szp1 ? Matrix{ComplexF64}(undef, Nvck_Szp1, Nvck_Szp1) : nothing
    bandSz0 = has_Sz0 ? Matrix{Float64}(undef, Nvck_Sz0, Nq) : nothing
    bandSzn1 = has_Szn1 ? Matrix{Float64}(undef, Nvck_Szn1, Nq) : nothing
    bandSzp1 = has_Szp1 ? Matrix{Float64}(undef, Nvck_Szp1, Nq) : nothing
    isallocated_bse = isallocated(bse)
    if !isallocated_bse
        start_time = time()
        allocate!(bse)
        allocate_elapsed = time() - start_time
        println("`allocate!` costs $(round(allocate_elapsed, digits=1))s, ")
        flush(stdout)
    end

    for (iq, q) in enumerate(qpoints)
        start_time = time()
        (q, isΓ) = _BSE_preprocess_Γ(q, false, bse.period)
        bandkq_up, bandkq_dn = _BSE_preprocess_q!(bse, q; Sz, isΓ)
        HSz0_view = has_Sz0 ? _BSE_Hamiltonian!(Val(0), HSz0, bse, q; isΓ) : nothing
        HSzn1_view = has_Szn1 ? _BSE_Hamiltonian!(Val(-1), HSzn1, bse) : nothing
        HSzp1_view = has_Szp1 ? _BSE_Hamiltonian!(Val(1), HSzp1, bse) : nothing
        Hamiltonian_time = time()
        has_Sz0 && (bandSz0[:, iq] = eigvals!(HSz0_view))
        has_Szn1 && (bandSzn1[:, iq] = eigvals!(HSzn1_view))
        has_Szp1 && (bandSzp1[:, iq] = eigvals!(HSzp1_view))
        eigen_time = time()
        Hamiltonian_elapsed = Hamiltonian_time - start_time
        eigen_elapsed = eigen_time - Hamiltonian_time
        println("$iq: $q; ",
            "Hamiltonian costs $(round(Hamiltonian_elapsed, digits=1))s, ",
            "eigen costs $(round(eigen_elapsed, digits=1))s.")
        flush(stdout)
    end
    !isallocated_bse && release!(bse)
    return (bandSz0=bandSz0, bandSzn1=bandSzn1, bandSzp1=bandSzp1)
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEspinaware_wiTDA, ::Val{false}, ::Val{:Periodic})
    return BAND_BSE(LinearAlgebra_BSEeigenStrategy, qpoints, bse, Val(false), Val(:Bloch))
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEspinaware_wiTDA, ::Val{false}, ::Val{:BlochPeriodic})
    return BAND_BSE(LinearAlgebra_BSEeigenStrategy, qpoints, bse, Val(false), Val(:Bloch))
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEspinaware_wiTDA, ::Val{true}, ::Val{:Bloch})

    Nq = length(qpoints)
    (Nvck_Sz0, Nvck_Szn1, Nvck_Szp1) = bse(Val(:dimension))
    has_Sz0 = !iszero(Nvck_Sz0)
    has_Szn1 = !iszero(Nvck_Szn1)
    has_Szp1 = !iszero(Nvck_Szp1)
    Sz = Int[]
    has_Sz0 && push!(Sz, 0)
    has_Szn1 && push!(Sz, -1)
    has_Szp1 && push!(Sz, 1)
    HSz0 = has_Sz0 ? Matrix{ComplexF64}(undef, Nvck_Sz0, Nvck_Sz0) : nothing
    HSzn1 = has_Szn1 ? Matrix{ComplexF64}(undef, Nvck_Szn1, Nvck_Szn1) : nothing
    HSzp1 = has_Szp1 ? Matrix{ComplexF64}(undef, Nvck_Szp1, Nvck_Szp1) : nothing
    bandSz0 = has_Sz0 ? Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}(undef, Nq) : nothing
    bandSzn1 = has_Szn1 ? Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}(undef, Nq) : nothing
    bandSzp1 = has_Szp1 ? Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}(undef, Nq) : nothing
    isallocated_bse = isallocated(bse)
    if !isallocated_bse
        start_time = time()
        allocate!(bse)
        allocate_elapsed = time() - start_time
        println("`allocate!` costs $(round(allocate_elapsed, digits=1))s, ")
        flush(stdout)
    end

    for (iq, q) in enumerate(qpoints)
        start_time = time()
        (q, isΓ) = _BSE_preprocess_Γ(q, false, bse.period)
        bandkq_up, bandkq_dn = _BSE_preprocess_q!(bse, q; Sz, isΓ)
        HSz0_view = has_Sz0 ? _BSE_Hamiltonian!(Val(0), HSz0, bse, q; isΓ) : nothing
        HSzn1_view = has_Szn1 ? _BSE_Hamiltonian!(Val(-1), HSzn1, bse) : nothing
        HSzp1_view = has_Szp1 ? _BSE_Hamiltonian!(Val(1), HSzp1, bse) : nothing
        Hamiltonian_time = time()
        has_Sz0 && (bandSz0[iq] = eigen!(HSz0_view))
        has_Szn1 && (bandSzn1[iq] = eigen!(HSzn1_view))
        has_Szp1 && (bandSzp1[iq] = eigen!(HSzp1_view))
        eigen_time = time()
        Hamiltonian_elapsed = Hamiltonian_time - start_time
        eigen_elapsed = eigen_time - Hamiltonian_time
        println("$iq: $q; ",
            "Hamiltonian costs $(round(Hamiltonian_elapsed, digits=1))s, ",
            "eigen costs $(round(eigen_elapsed, digits=1))s.")
        flush(stdout)
    end
    !isallocated_bse && release!(bse)
    return (bandSz0=bandSz0, bandSzn1=bandSzn1, bandSzp1=bandSzp1)
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEspinaware_wiTDA, ::Val{true}, ::Val{:Periodic})

    e_kR = [cispi(-2 * (k ⋅ R)) for R in bse.Rgrid, k in bse.kgrid]

    Nq = length(qpoints)
    (Nvck_Sz0, Nvck_Szn1, Nvck_Szp1) = bse(Val(:dimension))
    has_Sz0 = !iszero(Nvck_Sz0)
    has_Szn1 = !iszero(Nvck_Szn1)
    has_Szp1 = !iszero(Nvck_Szp1)
    Sz = Int[]
    has_Sz0 && push!(Sz, 0)
    has_Szn1 && push!(Sz, -1)
    has_Szp1 && push!(Sz, 1)
    HSz0 = has_Sz0 ? Matrix{ComplexF64}(undef, Nvck_Sz0, Nvck_Sz0) : nothing
    HSzn1 = has_Szn1 ? Matrix{ComplexF64}(undef, Nvck_Szn1, Nvck_Szn1) : nothing
    HSzp1 = has_Szp1 ? Matrix{ComplexF64}(undef, Nvck_Szp1, Nvck_Szp1) : nothing
    ubandSz0 = has_Sz0 ? Vector{BSE_uEigen}(undef, Nq) : nothing
    ubandSzn1 = has_Szn1 ? Vector{BSE_uEigen}(undef, Nq) : nothing
    ubandSzp1 = has_Szp1 ? Vector{BSE_uEigen}(undef, Nq) : nothing
    isallocated_bse = isallocated(bse)
    if !isallocated_bse
        start_time = time()
        allocate!(bse)
        allocate_elapsed = time() - start_time
        println("`allocate!` costs $(round(allocate_elapsed, digits=1))s, ")
        flush(stdout)
    end

    for (iq, q) in enumerate(qpoints)
        start_time = time()
        (q, isΓ) = _BSE_preprocess_Γ(q, false, bse.period)
        bandkq_up, bandkq_dn = _BSE_preprocess_q!(bse, q; Sz, isΓ)
        HSz0_view = has_Sz0 ? _BSE_Hamiltonian!(Val(0), HSz0, bse, q; isΓ) : nothing
        HSzn1_view = has_Szn1 ? _BSE_Hamiltonian!(Val(-1), HSzn1, bse) : nothing
        HSzp1_view = has_Szp1 ? _BSE_Hamiltonian!(Val(1), HSzp1, bse) : nothing
        Hamiltonian_time = time()
        has_Sz0 && (ubandSz0[iq] = _BSE_uEigen(bse, q, eigen!(HSz0_view), bandkq_up, bandkq_dn, e_kR, Val(0)))
        has_Szn1 && (ubandSzn1[iq] = _BSE_uEigen(bse, q, eigen!(HSzn1_view), bandkq_up, bandkq_dn, e_kR, Val(-1)))
        has_Szp1 && (ubandSzp1[iq] = _BSE_uEigen(bse, q, eigen!(HSzp1_view), bandkq_up, bandkq_dn, e_kR, Val(1)))
        eigen_time = time()
        Hamiltonian_elapsed = Hamiltonian_time - start_time
        eigen_elapsed = eigen_time - Hamiltonian_time
        println("$iq: $q; ",
            "Hamiltonian costs $(round(Hamiltonian_elapsed, digits=1))s, ",
            "eigen costs $(round(eigen_elapsed, digits=1))s.")
        flush(stdout)
    end
    !isallocated_bse && release!(bse)
    return (ubandSz0=ubandSz0, ubandSzn1=ubandSzn1, ubandSzp1=ubandSzp1)
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEspinaware_wiTDA, ::Val{true}, ::Val{:BlochPeriodic})

    e_kR = [cispi(-2 * (k ⋅ R)) for R in bse.Rgrid, k in bse.kgrid]

    Nq = length(qpoints)
    (Nvck_Sz0, Nvck_Szn1, Nvck_Szp1) = bse(Val(:dimension))
    has_Sz0 = !iszero(Nvck_Sz0)
    has_Szn1 = !iszero(Nvck_Szn1)
    has_Szp1 = !iszero(Nvck_Szp1)
    Sz = Int[]
    has_Sz0 && push!(Sz, 0)
    has_Szn1 && push!(Sz, -1)
    has_Szp1 && push!(Sz, 1)
    HSz0 = has_Sz0 ? Matrix{ComplexF64}(undef, Nvck_Sz0, Nvck_Sz0) : nothing
    HSzn1 = has_Szn1 ? Matrix{ComplexF64}(undef, Nvck_Szn1, Nvck_Szn1) : nothing
    HSzp1 = has_Szp1 ? Matrix{ComplexF64}(undef, Nvck_Szp1, Nvck_Szp1) : nothing
    bandSz0 = has_Sz0 ? Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}(undef, Nq) : nothing
    bandSzn1 = has_Szn1 ? Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}(undef, Nq) : nothing
    bandSzp1 = has_Szp1 ? Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}(undef, Nq) : nothing
    ubandSz0 = has_Sz0 ? Vector{BSE_uEigen}(undef, Nq) : nothing
    ubandSzn1 = has_Szn1 ? Vector{BSE_uEigen}(undef, Nq) : nothing
    ubandSzp1 = has_Szp1 ? Vector{BSE_uEigen}(undef, Nq) : nothing
    isallocated_bse = isallocated(bse)
    if !isallocated_bse
        start_time = time()
        allocate!(bse)
        allocate_elapsed = time() - start_time
        println("`allocate!` costs $(round(allocate_elapsed, digits=1))s, ")
        flush(stdout)
    end

    for (iq, q) in enumerate(qpoints)
        start_time = time()
        (q, isΓ) = _BSE_preprocess_Γ(q, false, bse.period)
        bandkq_up, bandkq_dn = _BSE_preprocess_q!(bse, q; Sz, isΓ)
        HSz0_view = has_Sz0 ? _BSE_Hamiltonian!(Val(0), HSz0, bse, q; isΓ) : nothing
        HSzn1_view = has_Szn1 ? _BSE_Hamiltonian!(Val(-1), HSzn1, bse) : nothing
        HSzp1_view = has_Szp1 ? _BSE_Hamiltonian!(Val(1), HSzp1, bse) : nothing
        Hamiltonian_time = time()
        has_Sz0 && (bandSz0[iq] = eigen!(HSz0_view))
        has_Szn1 && (bandSzn1[iq] = eigen!(HSzn1_view))
        has_Szp1 && (bandSzp1[iq] = eigen!(HSzp1_view))
        has_Sz0 && (ubandSz0[iq] = _BSE_uEigen(bse, q, bandSz0[iq], bandkq_up, bandkq_dn, e_kR, Val(0)))
        has_Szn1 && (ubandSzn1[iq] = _BSE_uEigen(bse, q, bandSzn1[iq], bandkq_up, bandkq_dn, e_kR, Val(-1)))
        has_Szp1 && (ubandSzp1[iq] = _BSE_uEigen(bse, q, bandSzp1[iq], bandkq_up, bandkq_dn, e_kR, Val(1)))
        eigen_time = time()
        Hamiltonian_elapsed = Hamiltonian_time - start_time
        eigen_elapsed = eigen_time - Hamiltonian_time
        println("$iq: $q; ",
            "Hamiltonian costs $(round(Hamiltonian_elapsed, digits=1))s, ",
            "eigen costs $(round(eigen_elapsed, digits=1))s.")
        flush(stdout)
    end
    !isallocated_bse && release!(bse)
    return (bandSz0=bandSz0, bandSzn1=bandSzn1, bandSzp1=bandSzp1,
        ubandSz0=ubandSz0, ubandSzn1=ubandSzn1, ubandSzp1=ubandSzp1)
end
