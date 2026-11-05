struct BSEspinblind_wiTDA{
    eleHT<:HermitianReciprocalHoppings,
    KᵈT<:PairInterAction,
    KˣT<:PairInterAction,
    bufferT<:_BSE_buffer
} <: AbstractBSE
    lattice::Lattice{Float64}
    atom_name::Vector{String}
    atom_location::Vector{ReducedCoordinates{Float64}}
    orb_name::Vector{String}
    orb_location::Vector{ReducedCoordinates{Float64}}
    period::Vec3{Bool}
    eleH::eleHT
    scissor::Float64
    Kᵈ::KᵈT
    Kˣ::KˣT
    channel::Int #the number of exciting channels. If 1, then output mixed; if 2, then output singlet and triplet.
    upindex::Vector{Int}
    dnindex::Vector{Int}
    kgrid::RedKgrid
    Rgrid::Vector{ReducedCoordinates{Int}}
    ndegen::Vector{Int}
    vck::TripleList{Int,Int,Int}
    ijR::TripleList{Int,Int,Int}
    bandk::Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}
    berryconnection::Matrix{ComplexF64} # [3, nvck]
    buffer::bufferT
end
mutable struct _BSEspinblind_wiTDA_buffer <: _BSE_buffer
    isallocated::Bool
    ΔE::Vector{Float64} # [nvck]
    Kᵈ_ij::Vector{NTuple{2,Int}} # [nij_W]
    ekq_pR::Matrix{ComplexF64} # [npR_W, nk]
    Kᵈ_pair::Vector{NTuple{2,Int}} # [np_W]
    cckk::Vector{NTuple{4,Int}} # [nc, nc, nk, nk]
    Mcckk::Array{ComplexF64,5} # [np_W, nc, nc, nk, nk]
    WMvvkk::Array{ComplexF64,5} # [np_W, nv, nv, nk, nk]
    Kˣ_ij::Vector{NTuple{2,Int}} # [nij_V]
    ek_pR::Matrix{ComplexF64} # [npR_V, nk]
    Kˣ_pair::Vector{NTuple{2,Int}} # [np_V]
    Mvck::Matrix{ComplexF64} # [np_V, nvck]
    Vq::Matrix{ComplexF64} # [np_V, np_V]
    VMvck::Matrix{ComplexF64} # [np_V, nvck]
end
function Base.empty(::Type{_BSEspinblind_wiTDA_buffer})
    return _BSEspinblind_wiTDA_buffer(false,
        Vector{Float64}(undef, 0),
        Vector{NTuple{2,Int}}(undef, 0),
        Matrix{ComplexF64}(undef, 0, 0),
        Vector{NTuple{2,Int}}(undef, 0),
        Vector{NTuple{4,Int}}(undef, 0),
        Array{ComplexF64,5}(undef, 0, 0, 0, 0, 0),
        Array{ComplexF64,5}(undef, 0, 0, 0, 0, 0),
        Vector{NTuple{2,Int}}(undef, 0),
        Matrix{ComplexF64}(undef, 0, 0),
        Vector{NTuple{2,Int}}(undef, 0),
        Matrix{ComplexF64}(undef, 0, 0),
        Matrix{ComplexF64}(undef, 0, 0),
        Matrix{ComplexF64}(undef, 0, 0))
end
function Base.empty!(buffer::_BSEspinblind_wiTDA_buffer)
    buffer.isallocated = false
    buffer.ΔE = Vector{Float64}(undef, 0)
    buffer.Kᵈ_ij = Vector{NTuple{2,Int}}(undef, 0)
    buffer.ekq_pR = Matrix{ComplexF64}(undef, 0, 0)
    buffer.Kᵈ_pair = Vector{NTuple{2,Int}}(undef, 0)
    buffer.cckk = Vector{NTuple{4,Int}}(undef, 0)
    buffer.Mcckk = Array{ComplexF64,5}(undef, 0, 0, 0, 0, 0)
    buffer.WMvvkk = Array{ComplexF64,5}(undef, 0, 0, 0, 0, 0)
    buffer.Kˣ_ij = Vector{NTuple{2,Int}}(undef, 0)
    buffer.ek_pR = Matrix{ComplexF64}(undef, 0, 0)
    buffer.Kˣ_pair = Vector{NTuple{2,Int}}(undef, 0)
    buffer.Mvck = Matrix{ComplexF64}(undef, 0, 0)
    buffer.Vq = Matrix{ComplexF64}(undef, 0, 0)
    buffer.VMvck = Matrix{ComplexF64}(undef, 0, 0)
    return buffer
end
WannierInterActionBase.numorb(bse::BSEspinblind_wiTDA) = length(bse.orb_location)
function Base.show(io::IO, bse::BSEspinblind_wiTDA)
    print(io, "$(count(bse.period)) dimensinal BSE model with $(numatom(bse)) atoms and $(numorb(bse)) orbitals.")
end
"""
    BSE_wiTDA(kgrid, TB::AbstractTightBindModel, Kᵈ::PairInterAction, Kˣ::PairInterAction; 
    v::AbstractVector{<:Integer}, c::AbstractVector{<:Integer}, scissor::Real=0, channel::Integer=2,
    upindex::AbstractVector{<:Integer}=Int[], dnindex::AbstractVector{<:Integer}=Int[])
"""
function BSE_wiTDA(kgrid::RedKgrid, TB::AbstractTightBindModel, Kᵈ::PairInterAction, Kˣ::PairInterAction;
    v::AbstractVector{<:Integer}, c::AbstractVector{<:Integer}, scissor::Real=0, channel::Integer=2,
    upindex::AbstractVector{<:Integer}=Int[], dnindex::AbstractVector{<:Integer}=Int[])

    @assert length(v) > 0 "Empty `v` is unaccepted!"
    @assert length(c) > 0 "Empty `c` is unaccepted!"
    nw = numorb(TB)
    @assert all(i->1 ≤ i ≤ nw, v) "all `v` must be in the range of [1, $nw]!"
    @assert all(i->1 ≤ i ≤ nw, c) "all `c` must be in the range of [1, $nw]!"
    @assert channel == 1 || channel == 2 "channel must be 1 or 2!"
    if channel == 1 # corresponding to spinor or single-spin.
        if Kᵈ.nw == Kˣ.nw
            if Kᵈ.nw * 2 == nw
                upindex, dnindex = _BSEspinblind_process_spinindex(upindex, dnindex, Kᵈ.nw)
            elseif Kᵈ.nw ≠ nw
                error("Incompatible dimensions between TB and interaction!")
            end
        elseif Kᵈ.nw ≠ Kˣ.nw
            if Kᵈ.nw * 2 == nw && Kˣ.nw == nw
                upindex, dnindex = _BSEspinblind_process_spinindex(upindex, dnindex, Kᵈ.nw)
            elseif Kᵈ.nw == nw && Kˣ.nw * 2 == nw
                upindex, dnindex = _BSEspinblind_process_spinindex(upindex, dnindex, Kˣ.nw)
            else
                error("Incompatible dimensions between TB and interaction!")
            end
        end
    elseif Kᵈ.nw ≠ nw || Kˣ.nw ≠ nw
        error("Incompatible dimensions between TB and interaction!")
    end

    # Rgrid
    Rgrid, ndegen = gridindex(Val(:WignerSeitz), TB.lattice, kgrid.kgrid_size)

    nk = length(kgrid)
    nR = length(Rgrid)

    v = collect(Int, v)
    c = collect(Int, c)
    vck = TripleList(v, c, collect(1:nk))
    ijR = TripleList(collect(1:nw), collect(1:nw), collect(1:nR))

    bandk = BAND(kgrid, TB.H; vector=true)
    phase_normalize!(:positive_maxabs, bandk)
    berryconnection = _BSE_berryconnection_vck(TB.H, kgrid, bandk, TB.lattice, TB.orb_location, vck)

    buffer = empty(_BSEspinblind_wiTDA_buffer)
    return BSEspinblind_wiTDA{typeof(TB.H),typeof(Kᵈ),typeof(Kˣ),typeof(buffer)}(
        TB.lattice, TB.atom_name, TB.atom_location, TB.orb_name, TB.orb_location, TB.period, TB.H, scissor,
        Kᵈ, Kˣ, channel, upindex, dnindex, kgrid, Rgrid, ndegen, vck, ijR,
        bandk, berryconnection, buffer)
end
function _BSEspinblind_process_spinindex(upindex, dnindex, nw)
    has_up = !isempty(upindex)
    has_dn = !isempty(dnindex)
    has_up || has_dn ||
        error("for spinor TB with spinless interaction, at least one of `upindex` or `dnindex` must be provided")
    if has_up
        @assert length(upindex) == nw "length of `upindex` must be `numorb(TB)/2`"
        @assert allunique(upindex) "`upindex` cannot contain duplicate elements"
        @assert all(i -> 1 ≤ i ≤ 2nw, upindex) "all entries of `upindex` must lie in [1, $(2nw)]"
    end
    if has_dn
        @assert length(dnindex) == nw "length of `dnindex` must be `numorb(TB)/2`"
        @assert allunique(dnindex) "`dnindex` cannot contain duplicate elements"
        @assert all(i -> 1 ≤ i ≤ 2nw, dnindex) "all entries of `dnindex` must lie in [1, $(2nw)]"
    end
    if has_up && has_dn
        @assert length(union(upindex, dnindex)) == 2nw "`upindex` and `dnindex` must be complementary"
        upindex = collect(Int, upindex)
        dnindex = collect(Int, dnindex)
    elseif has_up
        upindex = collect(Int, upindex)
        dnindex = setdiff(1:2nw, upindex)
    else
        dnindex = collect(Int, dnindex)
        upindex = setdiff(1:2nw, dnindex)
    end
    return upindex, dnindex
end
function allocate!(bse::BSEspinblind_wiTDA; force::Bool=false, storage::Symbol=:memory,
    Mmap_Mv::Bool=true, Mmap_Mc::Bool=false, Mmap_Mvck::Bool=false, Mmap_VMvck::Bool=false,
    Mmapfolder::AbstractString="./BSE.Mmap")

    if !force && isallocated(bse)
        return bse
    elseif force && isallocated(bse)
        release!(bse)
    else
        GC.gc()
    end

    vck = bse.vck
    Nvck = length(vck)
    v = values(vck)[1]
    c = values(vck)[2]
    kgrid = bse.kgrid
    nv = length(v)
    nc = length(c)
    nk = length(kgrid)
    Kᵈ = bse.Kᵈ
    Kˣ = bse.Kˣ

    np_W = Kᵈ.np
    np_V = Kˣ.np

    # report memory
    WMvvkk_size = format_memory(sizeof(ComplexF64) * np_W * nv * nv * nk * nk)
    Mcckk_size = format_memory(sizeof(ComplexF64) * np_W * nc * nc * nk * nk)
    Mvck_size = format_memory(sizeof(ComplexF64) * np_V * Nvck)
    VMvck_size = format_memory(sizeof(ComplexF64) * np_V * Nvck)
    println("WMvvkk: ComplexF64, dims=($np_W, $nv, $nv, $nk, $nk), total=$WMvvkk_size")
    println("Mcckk : ComplexF64, dims=($np_W, $nc, $nc, $nk, $nk), total=$Mcckk_size")
    println("Mvck  : ComplexF64, dims=($np_V, $Nvck), total=$Mvck_size")
    println("VMvck : ComplexF64, dims=($np_V, $Nvck), total=$VMvck_size")
    flush(stdout)

    if storage === :memory
        WMvvkk = Array{ComplexF64,5}(undef, np_W, nv, nv, nk, nk)
        Mcckk = Array{ComplexF64,5}(undef, np_W, nc, nc, nk, nk)
        Mvck = Matrix{ComplexF64}(undef, np_V, Nvck)
        VMvck = Matrix{ComplexF64}(undef, np_V, Nvck)
    elseif storage === :Mmap
        rm(Mmapfolder; recursive=true, force=true)
        mkpath(Mmapfolder)
        WMvvkk = Mmap_Mv && !iszero(np_W) ? Mmap.mmap(joinpath(Mmapfolder, "WMv"), Array{ComplexF64,5}, (np_W, nv, nv, nk, nk)) :
                 Array{ComplexF64,5}(undef, np_W, nv, nv, nk, nk)
        Mcckk = Mmap_Mc && !iszero(np_W) ? Mmap.mmap(joinpath(Mmapfolder, "Mc"), Array{ComplexF64,5}, (np_W, nc, nc, nk, nk)) :
                Array{ComplexF64,5}(undef, np_W, nc, nc, nk, nk)
        Mvck = Mmap_Mvck && !iszero(np_V) ? Mmap.mmap(joinpath(Mmapfolder, "Mvck"), Matrix{ComplexF64}, (np_V, Nvck)) :
               Matrix{ComplexF64}(undef, np_V, Nvck)
        VMvck = Mmap_VMvck && !iszero(np_V) ? Mmap.mmap(joinpath(Mmapfolder, "VMvck"), Matrix{ComplexF64}, (np_V, Nvck)) :
                Matrix{ComplexF64}(undef, np_V, Nvck)
    end

    Kᵈ_ij = unique(map(p->(p[1], p[2]), Kᵈ.pair))
    Kᵈ_ij_idx = Dict{NTuple{2,Int},Int}(ij=>i for (i, ij) in enumerate(Kᵈ_ij))
    Kᵈ_pair = map(p->(Kᵈ_ij_idx[(p[1], p[2])], p[3]), Kᵈ.pair)

    LinearAlgebra.BLAS.set_num_threads(1)
    nt = Threads.nthreads()

    # preprocess W
    kgrid_Γ = RedKgrid(MonkhorstPack(kgrid.kgrid_size))
    kgrid_minusmap = kgridmap(kgrid_Γ, kgrid, kgrid, -)
    iΓ = kgrid_minusmap[1, 1]
    Wk = Array{ComplexF64}(undef, np_W, np_W, nk)
    Threads.@threads for ik in Base.OneTo(nk)
        Wk_view = view(Wk, :, :, ik)
        Kᵈ(Wk_view, kgrid_Γ[ik], kgrid.kgrid_size; isΓ=(ik == iΓ))
        Wk_view ./= nk
    end
    # preprocess WMvvkk
    ek_pR = [cispi(-2 * (k ⋅ R)) for R in Kᵈ.pR, k in kgrid]
    vvkk = collect(Iterators.product(1:nv, 1:nv, 1:nk, 1:nk))
    vvkk_chunk = _split_tasks(vvkk, nt)
    # The internal variables of the multi-threaded loop are bound to the same-named variables outside the loop, 
    # thereby causing race conditions. 
    # Therefore, it is necessary to define the multi-threading loop as an independent function.
    _BSEspinblind_WMvvkk!(Val(Kᵈ.nw == numorb(bse)),
        WMvvkk, vvkk_chunk, bse.bandk, v, Kᵈ_ij, ek_pR, Kᵈ_pair, Wk, kgrid_minusmap,
        Kᵈ.nw, bse.upindex, bse.dnindex)

    LinearAlgebra.BLAS.set_num_threads(nt)

    Kˣ_ij = unique(map(p->(p[1], p[2]), Kˣ.pair))
    Kˣ_ij_idx = Dict{NTuple{2,Int},Int}(ij=>i for (i, ij) in enumerate(Kˣ_ij))
    Kˣ_pair = map(p->(Kˣ_ij_idx[(p[1], p[2])], p[3]), Kˣ.pair)

    buffer = bse.buffer
    buffer.isallocated = true
    buffer.ΔE = Vector{Float64}(undef, Nvck)
    buffer.Kᵈ_ij = Kᵈ_ij
    buffer.ekq_pR = Matrix{ComplexF64}(undef, length(Kᵈ.pR), nk)
    buffer.Kᵈ_pair = Kᵈ_pair
    buffer.cckk = vec(collect(Iterators.product(1:nc, 1:nc, 1:nk, 1:nk)))
    buffer.Mcckk = Mcckk
    buffer.WMvvkk = WMvvkk
    buffer.Kˣ_ij = Kˣ_ij
    buffer.ek_pR = [cispi(-2 * (k ⋅ R)) for R in Kˣ.pR, k in kgrid]
    buffer.Kˣ_pair = Kˣ_pair
    buffer.Mvck = Mvck
    buffer.Vq = Matrix{ComplexF64}(undef, np_V, np_V)
    buffer.VMvck = VMvck

    return bse
end
function release!(bse::BSEspinblind_wiTDA)
    isallocated(bse) || return bse
    empty!(bse.buffer)
    GC.gc()
    return bse
end
function _BSE_preprocess_q!(bse::BSEspinblind_wiTDA, q::ReducedCoordinates; isΓ::Bool)
    isallocated(bse) || allocate!(bse)
    if isΓ
        # q is direction
        q_value = ReducedCoordinates{Int}(0, 0, 0)
    else
        q_value = q
    end
    # bandkq
    bandk = bse.bandk
    if isΓ
        bandkq = bandk
    else
        bandkq = BAND(map(k -> k + q_value, bse.kgrid), bse.eleH; vector=true)
        phase_normalize!(:positive_maxabs, bandkq)
    end
    # ΔE
    vck = bse.vck
    scissor = bse.scissor
    ΔE = bse.buffer.ΔE
    for (ivck, (v, c, ik)) in enumerate(vck)
        ΔE[ivck] = bandkq[ik].values[c] - bandk[ik].values[v] + scissor
    end

    nt = Threads.nthreads()

    # Mcckk
    kgrid = bse.kgrid
    Kᵈ = bse.Kᵈ
    ekq_pR = bse.buffer.ekq_pR
    for ik in eachindex(kgrid)
        kq = kgrid[ik] + q_value
        for (ipR, R) in enumerate(Kᵈ.pR)
            ekq_pR[ipR, ik] = cispi(-2 * (kq ⋅ R))
        end
    end
    cckk_chunk = _split_tasks(bse.buffer.cckk, nt)
    # The internal variables of the multi-threaded loop are bound to the same-named variables outside the loop, 
    # thereby causing race conditions. 
    # Therefore, it is necessary to define the multi-threading loop as an independent function.
    _BSEspinblind_Mcckk!(Val(Kᵈ.nw == numorb(bse)),
        bse.buffer.Mcckk, cckk_chunk, bandkq, values(vck)[2], bse.buffer.Kᵈ_ij, ekq_pR, bse.buffer.Kᵈ_pair,
        Kᵈ.nw, bse.upindex, bse.dnindex)

    # Mvck
    Kˣ = bse.Kˣ
    Mvck = bse.buffer.Mvck
    ivck_chunk = _split_tasks(length(vck), nt)
    # The internal variables of the multi-threaded loop are bound to the same-named variables outside the loop, 
    # thereby causing race conditions. 
    # Therefore, it is necessary to define the multi-threading loop as an independent function.
    _BSEspinblind_Mvck!(Val(Kˣ.nw == numorb(bse)),
        Mvck, ivck_chunk, vck, bandk, bandkq, bse.buffer.Kˣ_ij, bse.buffer.ek_pR, bse.buffer.Kˣ_pair,
        Kˣ.nw, bse.upindex, bse.dnindex)
    # Vq
    Vq = bse.buffer.Vq
    VMvck = bse.buffer.VMvck
    Kˣ(Vq, q_value; isΓ)
    Vq ./= length(kgrid)
    mul!(VMvck, Vq, Mvck)

    return bandkq
end
function _BSEspinblind_WMvvkk!(::Val{true}, WMvvkk, vvkk_chunk, bandk, v, Kᵈ_ij, ek_pR, Kᵈ_pair, Wk, kgrid_minusmap,
    Kᵈ_nw, upindex, dnindex)
    nij_W = length(Kᵈ_ij)
    np_W = length(Kᵈ_pair)
    Threads.@threads for vvkks in vvkk_chunk
        buffer_nij_W = Vector{ComplexF64}(undef, nij_W)
        buffer_np_W = Vector{ComplexF64}(undef, np_W)
        @inbounds for (iv, iv′, ik, ik′) in vvkks
            ψ = view(bandk[ik].vectors, :, v[iv])
            ψ′ = view(bandk[ik′].vectors, :, v[iv′])
            @simd for ij in Base.OneTo(nij_W)
                (i, j) = Kᵈ_ij[ij]
                buffer_nij_W[ij] = ψ[i] * conj(ψ′[j])
            end
            ek_pR_k′ = view(ek_pR, :, ik′)
            @simd for ip in Base.OneTo(np_W)
                (ij, ipR) = Kᵈ_pair[ip]
                buffer_np_W[ip] = buffer_nij_W[ij] * ek_pR_k′[ipR]
            end
            W_kk′ = view(Wk, :, :, kgrid_minusmap[ik, ik′])
            WMvvkk_view = view(WMvvkk, :, iv, iv′, ik, ik′)
            mul!(WMvvkk_view, W_kk′, buffer_np_W)
        end
    end
end
function _BSEspinblind_WMvvkk!(::Val{false}, WMvvkk, vvkk_chunk, bandk, v, Kᵈ_ij, ek_pR, Kᵈ_pair, Wk, kgrid_minusmap,
    Kᵈ_nw, upindex, dnindex)
    nij_W = length(Kᵈ_ij)
    np_W = length(Kᵈ_pair)
    Threads.@threads for vvkks in vvkk_chunk
        buffer_nij_W = Vector{ComplexF64}(undef, nij_W)
        buffer_np_W = Vector{ComplexF64}(undef, np_W)
        ψ_u = Vector{ComplexF64}(undef, Kᵈ_nw)
        ψ_d = Vector{ComplexF64}(undef, Kᵈ_nw)
        ψ′_u = Vector{ComplexF64}(undef, Kᵈ_nw)
        ψ′_d = Vector{ComplexF64}(undef, Kᵈ_nw)
        @inbounds for (iv, iv′, ik, ik′) in vvkks
            vectors_ik, vv = bandk[ik].vectors, v[iv]
            ψ_u .= view(vectors_ik, upindex, vv)
            ψ_d .= view(vectors_ik, dnindex, vv)
            vectors_ik′, vv′ = bandk[ik′].vectors, v[iv′]
            ψ′_u .= view(vectors_ik′, upindex, vv′)
            ψ′_d .= view(vectors_ik′, dnindex, vv′)
            @simd for ij in Base.OneTo(nij_W)
                (i, j) = Kᵈ_ij[ij]
                buffer_nij_W[ij] = ψ_u[i] * conj(ψ′_u[j]) + ψ_d[i] * conj(ψ′_d[j])
            end
            ek_pR_k′ = view(ek_pR, :, ik′)
            @simd for ip in Base.OneTo(np_W)
                (ij, ipR) = Kᵈ_pair[ip]
                buffer_np_W[ip] = buffer_nij_W[ij] * ek_pR_k′[ipR]
            end
            W_kk′ = view(Wk, :, :, kgrid_minusmap[ik, ik′])
            WMvvkk_view = view(WMvvkk, :, iv, iv′, ik, ik′)
            mul!(WMvvkk_view, W_kk′, buffer_np_W)
        end
    end
end
function _BSEspinblind_Mcckk!(::Val{true}, Mcckk, cckk_chunk, bandkq, c, Kᵈ_ij, ekq_pR, Kᵈ_pair,
    Kᵈ_nw, upindex, dnindex)
    nij_W = length(Kᵈ_ij)
    np_W = length(Kᵈ_pair)
    Threads.@threads for cckks in cckk_chunk
        buffer_nij_W = Vector{ComplexF64}(undef, nij_W)
        @inbounds for (ic, ic′, ik, ik′) in cckks
            ψ = view(bandkq[ik].vectors, :, c[ic])
            ψ′ = view(bandkq[ik′].vectors, :, c[ic′])
            @simd for ij in Base.OneTo(nij_W)
                (i, j) = Kᵈ_ij[ij]
                buffer_nij_W[ij] = ψ[i] * conj(ψ′[j])
            end
            Mcckk_view = view(Mcckk, :, ic, ic′, ik, ik′)
            ekq_pR_k′ = view(ekq_pR, :, ik′)
            @simd for ip in Base.OneTo(np_W)
                (ij, ipR) = Kᵈ_pair[ip]
                Mcckk_view[ip] = buffer_nij_W[ij] * ekq_pR_k′[ipR]
            end
        end
    end
end
function _BSEspinblind_Mcckk!(::Val{false}, Mcckk, cckk_chunk, bandkq, c, Kᵈ_ij, ekq_pR, Kᵈ_pair,
    Kᵈ_nw, upindex, dnindex)
    nij_W = length(Kᵈ_ij)
    np_W = length(Kᵈ_pair)
    Threads.@threads for cckks in cckk_chunk
        buffer_nij_W = Vector{ComplexF64}(undef, nij_W)
        ψ_u = Vector{ComplexF64}(undef, Kᵈ_nw)
        ψ_d = Vector{ComplexF64}(undef, Kᵈ_nw)
        ψ′_u = Vector{ComplexF64}(undef, Kᵈ_nw)
        ψ′_d = Vector{ComplexF64}(undef, Kᵈ_nw)
        @inbounds for (ic, ic′, ik, ik′) in cckks
            vectors_ik, cc = bandkq[ik].vectors, c[ic]
            ψ_u .= view(vectors_ik, upindex, cc)
            ψ_d .= view(vectors_ik, dnindex, cc)
            vectors_ik′, cc′ = bandkq[ik′].vectors, c[ic′]
            ψ′_u .= view(vectors_ik′, upindex, cc′)
            ψ′_d .= view(vectors_ik′, dnindex, cc′)
            @simd for ij in Base.OneTo(nij_W)
                (i, j) = Kᵈ_ij[ij]
                buffer_nij_W[ij] = ψ_u[i] * conj(ψ′_u[j]) + ψ_d[i] * conj(ψ′_d[j])
            end
            Mcckk_view = view(Mcckk, :, ic, ic′, ik, ik′)
            ekq_pR_k′ = view(ekq_pR, :, ik′)
            @simd for ip in Base.OneTo(np_W)
                (ij, ipR) = Kᵈ_pair[ip]
                Mcckk_view[ip] = buffer_nij_W[ij] * ekq_pR_k′[ipR]
            end
        end
    end
end
function _BSEspinblind_Mvck!(::Val{true}, Mvck, ivck_chunk, vck, bandk, bandkq, Kˣ_ij, ek_pR, Kˣ_pair,
    Kˣ_nw, upindex, dnindex)
    nij_V = length(Kˣ_ij)
    np_V = length(Kˣ_pair)
    Threads.@threads for ivcks in ivck_chunk
        buffer_nij_V = Vector{ComplexF64}(undef, nij_V)
        @inbounds for ivck in ivcks
            (v, c, ik) = vck[ivck]
            ψc = view(bandkq[ik].vectors, :, c)
            ψv = view(bandk[ik].vectors, :, v)
            @simd for ij in Base.OneTo(nij_V)
                (i, j) = Kˣ_ij[ij]
                buffer_nij_V[ij] = ψc[i] * conj(ψv[j])
            end
            Mvck_view = view(Mvck, :, ivck)
            ek_pR_k = view(ek_pR, :, ik)
            @simd for ip in Base.OneTo(np_V)
                (ij, ipR) = Kˣ_pair[ip]
                Mvck_view[ip] = buffer_nij_V[ij] * ek_pR_k[ipR]
            end
        end
    end
end
function _BSEspinblind_Mvck!(::Val{false}, Mvck, ivck_chunk, vck, bandk, bandkq, Kˣ_ij, ek_pR, Kˣ_pair,
    Kˣ_nw, upindex, dnindex)
    nij_V = length(Kˣ_ij)
    np_V = length(Kˣ_pair)
    Threads.@threads for ivcks in ivck_chunk
        buffer_nij_V = Vector{ComplexF64}(undef, nij_V)
        ψcu = Vector{ComplexF64}(undef, Kˣ_nw)
        ψcd = Vector{ComplexF64}(undef, Kˣ_nw)
        ψvu = Vector{ComplexF64}(undef, Kˣ_nw)
        ψvd = Vector{ComplexF64}(undef, Kˣ_nw)
        @inbounds for ivck in ivcks
            (v, c, ik) = vck[ivck]
            vectors_ikq = bandkq[ik].vectors
            ψcu .= view(vectors_ikq, upindex, c)
            ψcd .= view(vectors_ikq, dnindex, c)
            vectors_ik = bandk[ik].vectors
            ψvu .= view(vectors_ik, upindex, v)
            ψvd .= view(vectors_ik, dnindex, v)
            @simd for ij in Base.OneTo(nij_V)
                (i, j) = Kˣ_ij[ij]
                buffer_nij_V[ij] = ψcu[i] * conj(ψvu[j]) + ψcd[i] * conj(ψvd[j])
            end
            Mvck_view = view(Mvck, :, ivck)
            ek_pR_k = view(ek_pR, :, ik)
            @simd for ip in Base.OneTo(np_V)
                (ij, ipR) = Kˣ_pair[ip]
                Mvck_view[ip] = buffer_nij_V[ij] * ek_pR_k[ipR]
            end
        end
    end
end
(bse::BSEspinblind_wiTDA)(::Val{:dimension}) = length(bse.vck)
# isΓ表示是否精确处理Γ点，这时输入的q表示方向，注意是分数坐标。
function (bse::BSEspinblind_wiTDA)(q::ReducedCoordinates; isΓ::Bool=false)
    (q, isΓ) = _BSE_preprocess_Γ(q, isΓ, bse.period)
    bandkq = _BSE_preprocess_q!(bse, q; isΓ)
    Nvck = length(bse.vck)
    if bse.channel == 1
        H = Matrix{ComplexF64}(undef, Nvck, Nvck)
        H = _BSE_Hamiltonian!(Val(:mixed), H, bse, q; isΓ)
        return (bandkq=bandkq, H=H)
    elseif bse.channel == 2
        Ht = Matrix{ComplexF64}(undef, Nvck, Nvck)
        Hs = Matrix{ComplexF64}(undef, Nvck, Nvck)
        Ht, Hs = _BSE_Hamiltonian!(Val(:separate), Ht, Hs, bse, q; isΓ)
        return (bandkq=bandkq, Ht=Ht, Hs=Hs)
    end
end
function (bse::BSEspinblind_wiTDA)(::Val{:triplet}, q::ReducedCoordinates; isΓ::Bool=false)
    @assert bse.channel == 2 "For channel=1, please use `bse(q; isqgrid, isΓ)`!"
    Nvck = length(bse.vck)
    Ht = Matrix{ComplexF64}(undef, Nvck, Nvck)
    (q, isΓ) = _BSE_preprocess_Γ(q, isΓ, bse.period)
    bandkq = _BSE_preprocess_q!(bse, q; isΓ)
    Ht = _BSE_Hamiltonian!(Val(:triplet), Ht, bse, q; isΓ)
    return (bandkq=bandkq, Ht=Ht)
end
function (bse::BSEspinblind_wiTDA)(::Val{:singlet}, q::ReducedCoordinates; isΓ::Bool=false)
    @assert bse.channel == 2 "For channel=1, please use `bse(q; isqgrid, isΓ)`!"
    Nvck = length(bse.vck)
    Hs = Matrix{ComplexF64}(undef, Nvck, Nvck)
    (q, isΓ) = _BSE_preprocess_Γ(q, isΓ, bse.period)
    bandkq = _BSE_preprocess_q!(bse, q; isΓ)
    Hs = _BSE_Hamiltonian!(Val(:singlet), Hs, bse, q; isΓ)
    return (bandkq=bandkq, Hs=Hs)
end
function (bse::BSEspinblind_wiTDA)(H::AbstractMatrix, q::ReducedCoordinates; isΓ::Bool=false)
    @assert bse.channel == 1 "For channel=2, please use `bse((Ht=Ht, Hs=Hs), q; isqgrid, isΓ)`!"
    Nvck = length(bse.vck)
    size(H) == (Nvck, Nvck) ||
        throw(DimensionMismatch("Hamiltonian must be a $Nvck×$Nvck matrix, got size $(size(H))"))
    (q, isΓ) = _BSE_preprocess_Γ(q, isΓ, bse.period)
    bandkq = _BSE_preprocess_q!(bse, q; isΓ)
    H = _BSE_Hamiltonian!(Val(:mixed), H, bse, q; isΓ)
    return (bandkq=bandkq, H=H)
end
function (bse::BSEspinblind_wiTDA)(Hlist::NamedTuple, q::ReducedCoordinates; isΓ::Bool=false)
    @assert issubset(keys(Hlist), (:H, :Ht, :Hs)) || "Hlist only can contain :H, :Ht and :Hs"
    Nvck = length(bse.vck)
    if bse.channel == 1
        H = get(Hlist, :H, nothing)
        has_H = H isa AbstractMatrix{<:Number}
        if has_H
            size(H) == (Nvck, Nvck) ||
                throw(DimensionMismatch("Hamiltonian must be a $Nvck×$Nvck matrix, got size $(size(H))"))
            (q, isΓ) = _BSE_preprocess_Γ(q, isΓ, bse.period)
            bandkq = _BSE_preprocess_q!(bse, q; isΓ)
            H = _BSE_Hamiltonian!(Val(:mixed), H, bse, q; isΓ)
        else
            bandkq = nothing
            H = nothing
        end
        return (bandkq=bandkq, H=H)
    else
        Ht = get(Hlist, :Ht, nothing)
        Hs = get(Hlist, :Hs, nothing)
        has_Ht = Ht isa AbstractMatrix{<:Number}
        has_Hs = Hs isa AbstractMatrix{<:Number}
        if has_Ht
            size(Ht) == (Nvck, Nvck) ||
                throw(DimensionMismatch("Htriplet must be a $Nvck×$Nvck matrix, got size $(size(Ht))"))
        end
        if has_Hs
            size(Hs) == (Nvck, Nvck) ||
                throw(DimensionMismatch("Hsinglet must be a $Nvck×$Nvck matrix, got size $(size(Hs))"))
        end
        if !has_Ht && !has_Hs
            return (bandkq=nothing, Ht=nothing, Hs=nothing)
        end
        (q, isΓ) = _BSE_preprocess_Γ(q, isΓ, bse.period)
        bandkq = _BSE_preprocess_q!(bse, q; isΓ)
        if has_Ht && has_Hs
            Ht, Hs = _BSE_Hamiltonian!(Val(:separate), Ht, Hs, bse, q; isΓ)
        elseif has_Ht
            Ht = _BSE_Hamiltonian!(Val(:triplet), Ht, bse, q; isΓ)
            Hs = nothing
        else
            Hs = _BSE_Hamiltonian!(Val(:singlet), Hs, bse, q; isΓ)
            Ht = nothing
        end
        return (bandkq=bandkq, Ht=Ht, Hs=Hs)
    end
end

function _BSE_Hamiltonian!(::Val{:mixed}, H::AbstractMatrix, bse::BSEspinblind_wiTDA, q; isΓ)
    vck = bse.vck
    ΔE = bse.buffer.ΔE
    Mcckk = bse.buffer.Mcckk
    WMvvkk = bse.buffer.WMvvkk
    Mvck = bse.buffer.Mvck
    VMvck = bse.buffer.VMvck
    LinearAlgebra.BLAS.set_num_threads(1)
    Threads.@threads for (ivck, ivck′) in UpperTriangularIndices(length(vck))
        (iv, ic, ik) = triple_index(vck, ivck)
        (iv′, ic′, ik′) = triple_index(vck, ivck′)
        Mc = view(Mcckk, :, ic, ic′, ik, ik′)
        WMv = view(WMvvkk, :, iv, iv′, ik, ik′)
        Kᵈ = - Mc ⋅ WMv
        M = view(Mvck, :, ivck)
        VM′ = view(VMvck, :, ivck′)
        Kˣ = M ⋅ VM′
        if ivck == ivck′
            H[ivck, ivck] = ΔE[ivck] + Kᵈ + Kˣ
        else
            H[ivck, ivck′] = Kᵈ + Kˣ
        end
    end
    # infinite term in Kˣ
    if isΓ && count(bse.period) == 3
        _BSE_KˣΓ!(Val(1), H, bse, q)
    end
    LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
    return Hermitian(H, :U)
end
function _BSE_Hamiltonian!(::Val{:triplet}, Ht::AbstractMatrix, bse::BSEspinblind_wiTDA, q; isΓ)
    vck = bse.vck
    ΔE = bse.buffer.ΔE
    Mcckk = bse.buffer.Mcckk
    WMvvkk = bse.buffer.WMvvkk
    LinearAlgebra.BLAS.set_num_threads(1)
    Threads.@threads for (ivck, ivck′) in UpperTriangularIndices(length(vck))
        (iv, ic, ik) = triple_index(vck, ivck)
        (iv′, ic′, ik′) = triple_index(vck, ivck′)
        Mc = view(Mcckk, :, ic, ic′, ik, ik′)
        WMv = view(WMvvkk, :, iv, iv′, ik, ik′)
        Kᵈ = - Mc ⋅ WMv
        if ivck == ivck′
            Ht[ivck, ivck] = ΔE[ivck] + Kᵈ
        else
            Ht[ivck, ivck′] = Kᵈ
        end
    end
    LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
    return Hermitian(Ht, :U)
end
function _BSE_Hamiltonian!(::Val{:singlet}, Hs::AbstractMatrix, bse::BSEspinblind_wiTDA, q; isΓ)
    vck = bse.vck
    ΔE = bse.buffer.ΔE
    Mcckk = bse.buffer.Mcckk
    WMvvkk = bse.buffer.WMvvkk
    Mvck = bse.buffer.Mvck
    VMvck = bse.buffer.VMvck
    LinearAlgebra.BLAS.set_num_threads(1)
    Threads.@threads for (ivck, ivck′) in UpperTriangularIndices(length(vck))
        (iv, ic, ik) = triple_index(vck, ivck)
        (iv′, ic′, ik′) = triple_index(vck, ivck′)
        Mc = view(Mcckk, :, ic, ic′, ik, ik′)
        WMv = view(WMvvkk, :, iv, iv′, ik, ik′)
        Kᵈ = - Mc ⋅ WMv
        M = view(Mvck, :, ivck)
        VM′ = view(VMvck, :, ivck′)
        Kˣ = M ⋅ VM′
        if ivck == ivck′
            Hs[ivck, ivck] = ΔE[ivck] + Kᵈ + 2 * Kˣ
        else
            Hs[ivck, ivck′] = Kᵈ + 2 * Kˣ
        end
    end
    # infinite term in Kˣ
    if isΓ && count(bse.period) == 3
        _BSE_KˣΓ!(Val(2), Hs, bse, q)
    end
    LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
    return Hermitian(Hs, :U)
end
function _BSE_Hamiltonian!(::Val{:separate}, Ht::AbstractMatrix, Hs::AbstractMatrix, bse::BSEspinblind_wiTDA, q; isΓ)
    vck = bse.vck
    ΔE = bse.buffer.ΔE
    Mcckk = bse.buffer.Mcckk
    WMvvkk = bse.buffer.WMvvkk
    Mvck = bse.buffer.Mvck
    VMvck = bse.buffer.VMvck
    LinearAlgebra.BLAS.set_num_threads(1)
    Threads.@threads for (ivck, ivck′) in UpperTriangularIndices(length(vck))
        (iv, ic, ik) = triple_index(vck, ivck)
        (iv′, ic′, ik′) = triple_index(vck, ivck′)
        Mc = view(Mcckk, :, ic, ic′, ik, ik′)
        WMv = view(WMvvkk, :, iv, iv′, ik, ik′)
        Kᵈ = - Mc ⋅ WMv
        M = view(Mvck, :, ivck)
        VM′ = view(VMvck, :, ivck′)
        Kˣ = M ⋅ VM′
        if ivck == ivck′
            Ht[ivck, ivck] = ΔE[ivck] + Kᵈ
            Hs[ivck, ivck] = ΔE[ivck] + Kᵈ + 2 * Kˣ
        else
            Ht[ivck, ivck′] = Kᵈ
            Hs[ivck, ivck′] = Kᵈ + 2 * Kˣ
        end
    end
    # infinite term in Kˣ
    if isΓ && count(bse.period) == 3
        _BSE_KˣΓ!(Val(2), Hs, bse, q)
    end
    LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
    return Hermitian(Ht, :U), Hermitian(Hs, :U)
end
function _BSE_KˣΓ!(::Val{1}, H::AbstractMatrix, bse::BSEspinblind_wiTDA, q)
    Ω = abs(det(parent(bse.lattice)))
    CoulombScaledivnk = CoulombScale * 4π / Ω / length(bse.kgrid)
    rlattice = reciprocal(bse.lattice)
    qcar = rlattice * q
    qnorm = norm(qcar)
    q̂x = qcar[1] / qnorm
    q̂y = qcar[2] / qnorm
    q̂z = qcar[3] / qnorm
    A = bse.berryconnection
    vck = bse.vck
    Threads.@threads for (ivck, ivck′) in UpperTriangularIndices(length(vck))
        H[ivck, ivck′] += conj(A[1, ivck] * q̂x + A[2, ivck] * q̂y + A[3, ivck] * q̂z) *
            (A[1, ivck′] * q̂x + A[2, ivck′] * q̂y + A[3, ivck′] * q̂z) * CoulombScaledivnk
    end
    return H
end
function _BSE_KˣΓ!(::Val{2}, Hs::AbstractMatrix, bse::BSEspinblind_wiTDA, q)
    Ω = abs(det(parent(bse.lattice)))
    CoulombScaledivnk = CoulombScale * 4π / Ω / length(bse.kgrid)
    rlattice = reciprocal(bse.lattice)
    qcar = rlattice * q
    qnorm = norm(qcar)
    q̂x = qcar[1] / qnorm
    q̂y = qcar[2] / qnorm
    q̂z = qcar[3] / qnorm
    A = bse.berryconnection
    vck = bse.vck
    Threads.@threads for (ivck, ivck′) in UpperTriangularIndices(length(vck))
        Hs[ivck, ivck′] += 2 * conj(A[1, ivck] * q̂x + A[2, ivck] * q̂y + A[3, ivck] * q̂z) *
            (A[1, ivck′] * q̂x + A[2, ivck′] * q̂y + A[3, ivck′] * q̂z) * CoulombScaledivnk
    end
    return Hs
end
