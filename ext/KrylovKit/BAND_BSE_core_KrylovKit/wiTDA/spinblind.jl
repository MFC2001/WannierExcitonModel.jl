function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEspinblind_wiTDA, ::Val{false}, ::Val{:Bloch})

    Nq = length(qpoints)
    Nvck = length(bse.vck)
    if bse.channel == 1
        Hlist = (H=Matrix{ComplexF64}(undef, Nvck, Nvck),)
        band=(band=Matrix{Float64}(undef, eigsolveconfig[:howmany], Nq),)
    else
        Hlist = (Ht=Matrix{ComplexF64}(undef, Nvck, Nvck), Hs=Matrix{ComplexF64}(undef, Nvck, Nvck))
        bandt = Matrix{Float64}(undef, eigsolveconfig[:howmany], Nq)
        bands = Matrix{Float64}(undef, eigsolveconfig[:howmany], Nq)
        band=(bandt=bandt, bands=bands)
    end
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
        result = bse(Hlist, q; isΓ=false)
        Hamiltonian_time = time()
        if bse.channel == 1
            band.band[:, iq] = _eigsolve_Hmat(result.H).values
        else
            band.bandt[:, iq] = _eigsolve_Hmat(result.Ht).values
            band.bands[:, iq] = _eigsolve_Hmat(result.Hs).values
        end
        eigen_time = time()
        Hamiltonian_elapsed = Hamiltonian_time - start_time
        eigen_elapsed = eigen_time - Hamiltonian_time
        println("$iq: $q; ",
            "Hamiltonian costs $(round(Hamiltonian_elapsed, digits=1))s, ",
            "eigen costs $(round(eigen_elapsed, digits=1))s.")
        flush(stdout)
    end
    !isallocated_bse && release!(bse)
    return band
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEspinblind_wiTDA, ::Val{false}, ::Val{:Periodic})
    return BAND_BSE(KrylovKit_BSEeigenStrategy, qpoints, bse, Val(false), Val(:Bloch))
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEspinblind_wiTDA, ::Val{false}, ::Val{:BlochPeriodic})
    return BAND_BSE(KrylovKit_BSEeigenStrategy, qpoints, bse, Val(false), Val(:Bloch))
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEspinblind_wiTDA, ::Val{true}, ::Val{:Bloch})

    Nq = length(qpoints)
    Nvck = length(bse.vck)
    if bse.channel == 1
        Hlist = (H=Matrix{ComplexF64}(undef, Nvck, Nvck),)
        band=(band=Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}(undef, Nq),)
    else
        Hlist = (Ht=Matrix{ComplexF64}(undef, Nvck, Nvck), Hs=Matrix{ComplexF64}(undef, Nvck, Nvck))
        bandt = Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}(undef, Nq)
        bands = Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}(undef, Nq)
        band=(bandt=bandt, bands=bands)
    end
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
        result = bse(Hlist, q; isΓ=false)
        Hamiltonian_time = time()
        if bse.channel == 1
            band.band[iq] = _eigsolve_Hmat(result.H)
        else
            band.bandt[iq] = _eigsolve_Hmat(result.Ht)
            band.bands[iq] = _eigsolve_Hmat(result.Hs)
        end
        eigen_time = time()
        Hamiltonian_elapsed = Hamiltonian_time - start_time
        eigen_elapsed = eigen_time - Hamiltonian_time
        println("$iq: $q; ",
            "Hamiltonian costs $(round(Hamiltonian_elapsed, digits=1))s, ",
            "eigen costs $(round(eigen_elapsed, digits=1))s.")
        flush(stdout)
    end
    !isallocated_bse && release!(bse)
    return band
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEspinblind_wiTDA, ::Val{true}, ::Val{:Periodic})

    Nq = length(qpoints)
    Nvck = length(bse.vck)
    if bse.channel == 1
        Hlist = (H=Matrix{ComplexF64}(undef, Nvck, Nvck),)
        band=(uband=Vector{BSEspinblind_wiTDA_uEigen}(undef, Nq),)
    else
        Hlist = (Ht=Matrix{ComplexF64}(undef, Nvck, Nvck), Hs=Matrix{ComplexF64}(undef, Nvck, Nvck))
        ubandt = Vector{BSEspinblind_wiTDA_uEigen}(undef, Nq)
        ubands = Vector{BSEspinblind_wiTDA_uEigen}(undef, Nq)
        band=(ubandt=ubandt, ubands=ubands)
    end
    e_kR = [cispi(-2 * (k ⋅ R)) for R in bse.Rgrid, k in bse.kgrid]
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
        result = bse(Hlist, q; isΓ=false)
        Hamiltonian_time = time()
        if bse.channel == 1
            band.uband[iq] = _BSE_uEigen(bse, q, _eigsolve_Hmat(result.H), result.bandkq, e_kR)
        else
            band.ubandt[iq] = _BSE_uEigen(bse, q, _eigsolve_Hmat(result.Ht), result.bandkq, e_kR)
            band.ubands[iq] = _BSE_uEigen(bse, q, _eigsolve_Hmat(result.Hs), result.bandkq, e_kR)
        end
        eigen_time = time()
        Hamiltonian_elapsed = Hamiltonian_time - start_time
        eigen_elapsed = eigen_time - Hamiltonian_time
        println("$iq: $q; ",
            "Hamiltonian costs $(round(Hamiltonian_elapsed, digits=1))s, ",
            "eigen costs $(round(eigen_elapsed, digits=1))s.")
        flush(stdout)
    end
    !isallocated_bse && release!(bse)
    return band
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEspinblind_wiTDA, ::Val{true}, ::Val{:BlochPeriodic})

    Nq = length(qpoints)
    Nvck = length(bse.vck)
    if bse.channel == 1
        Hlist = (H=Matrix{ComplexF64}(undef, Nvck, Nvck),)
        band=(band=Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}(undef, Nq),
            uband=Vector{BSEspinblind_wiTDA_uEigen}(undef, Nq))
    else
        Hlist = (Ht=Matrix{ComplexF64}(undef, Nvck, Nvck), Hs=Matrix{ComplexF64}(undef, Nvck, Nvck))
        bandt = Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}(undef, Nq)
        bands = Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}(undef, Nq)
        ubandt = Vector{BSEspinblind_wiTDA_uEigen}(undef, Nq)
        ubands = Vector{BSEspinblind_wiTDA_uEigen}(undef, Nq)
        band=(bandt=bandt, bands=bands, ubandt=ubandt, ubands=ubands)
    end
    e_kR = [cispi(-2 * (k ⋅ R)) for R in bse.Rgrid, k in bse.kgrid]
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
        result = bse(Hlist, q; isΓ=false)
        Hamiltonian_time = time()
        if bse.channel == 1
            band.band[iq] = _eigsolve_Hmat(result.H)
            band.uband[iq] = _BSE_uEigen(bse, q, band.band[iq], result.bandkq, e_kR)
        else
            band.bandt[iq] = _eigsolve_Hmat(result.Ht)
            band.bands[iq] = _eigsolve_Hmat(result.Hs)
            band.ubandt[iq] = _BSE_uEigen(bse, q, band.bandt[iq], result.bandkq, e_kR)
            band.ubands[iq] = _BSE_uEigen(bse, q, band.bands[iq], result.bandkq, e_kR)
        end
        eigen_time = time()
        Hamiltonian_elapsed = Hamiltonian_time - start_time
        eigen_elapsed = eigen_time - Hamiltonian_time
        println("$iq: $q; ",
            "Hamiltonian costs $(round(Hamiltonian_elapsed, digits=1))s, ",
            "eigen costs $(round(eigen_elapsed, digits=1))s.")
        flush(stdout)
    end
    !isallocated_bse && release!(bse)
    return band
end
