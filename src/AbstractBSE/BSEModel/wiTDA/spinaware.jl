struct BSEspinaware_wiTDA{
    eleHT_up<:HermitianReciprocalHoppings,
    eleHT_dn<:HermitianReciprocalHoppings,
    KᵈT<:PairInterAction,
    KˣT<:PairInterAction,
    bufferT<:_BSE_buffer
} <: AbstractBSE
    lattice::Lattice{Float64}
    atom_name::Vector{String}
    atom_location::Vector{ReducedCoordinates{Float64}}
    orb_name_up::Vector{String}
    orb_location_up::Vector{ReducedCoordinates{Float64}}
    orb_name_dn::Vector{String}
    orb_location_dn::Vector{ReducedCoordinates{Float64}}
    period::Vec3{Bool}
    eleH_up::eleHT_up
    eleH_dn::eleHT_dn
    scissor_c_up::Float64
    scissor_v_up::Float64
    scissor_c_dn::Float64
    scissor_v_dn::Float64
    Kᵈ::KᵈT
    Kˣ::KˣT
    kgrid::RedKgrid
    Rgrid::Vector{ReducedCoordinates{Int}}
    ndegen::Vector{Int}
    vck_uu::TripleList{Int,Int,Int}
    vck_dd::TripleList{Int,Int,Int}
    vck_ud::TripleList{Int,Int,Int}
    vck_du::TripleList{Int,Int,Int}
    ijR_uu::TripleList{Int,Int,Int}
    ijR_dd::TripleList{Int,Int,Int}
    ijR_ud::TripleList{Int,Int,Int}
    ijR_du::TripleList{Int,Int,Int}
    bandk_up::Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}
    bandk_dn::Vector{Eigen{ComplexF64,Float64,Matrix{ComplexF64},Vector{Float64}}}
    Sz::Set{Int} # 0, -1, 1
    berryconnection_up::Matrix{ComplexF64}
    berryconnection_dn::Matrix{ComplexF64}
    buffer::bufferT
end
mutable struct _BSEspinaware_wiTDA_Sz0_buffer <: _BSE_buffer
    ΔE_uu::Vector{Float64} # [nvck_uu]
    ΔE_dd::Vector{Float64} # [nvck_dd]
    WMvvkk_uu::Array{ComplexF64,5} # [np_W_up, nv_up, nv_up, nk, nk]
    WMvvkk_dd::Array{ComplexF64,5} # [np_W_dn, nv_dn, nv_dn, nk, nk]
    Kˣ_ij_up::Vector{NTuple{2,Int}} # [nij_V_up]
    Kˣ_ij_dn::Vector{NTuple{2,Int}} # [nij_V_dn]
    ek_pR::Matrix{ComplexF64} # [npR_V, nk]
    Kˣ_pair_up::Vector{NTuple{2,Int}} # [np_V_up]
    Kˣ_pair_dn::Vector{NTuple{2,Int}} # [np_V_dn]
    Mvck_uu::Matrix{ComplexF64} # [np_V_up, nvck_uu]
    Mvck_dd::Matrix{ComplexF64} # [np_V_dn, nvck_dd]
    Kˣ_pair_upindex::Vector{Int} # [np_V_up]
    Kˣ_pair_dnindex::Vector{Int} # [np_V_dn]
    Vq::Matrix{ComplexF64} # [np_V, np_V]
    Vq_uu::Matrix{ComplexF64} # [np_V_up, np_V_up]
    Vq_dd::Matrix{ComplexF64} # [np_V_dn, np_V_dn]
    Vq_ud::Matrix{ComplexF64} # [np_V_up, np_V_dn]
    VMvck_uu::Matrix{ComplexF64} # [np_V_up, nvck_uu]
    VMvck_dd::Matrix{ComplexF64} # [np_V_dn, nvck_dd]
    VMvck_ud::Matrix{ComplexF64} # [np_V_up, nvck_dd]
end
function Base.empty(::Type{_BSEspinaware_wiTDA_Sz0_buffer})
    return _BSEspinaware_wiTDA_Sz0_buffer(
        Vector{Float64}(undef, 0),
        Vector{Float64}(undef, 0),
        Array{ComplexF64,5}(undef, 0, 0, 0, 0, 0),
        Array{ComplexF64,5}(undef, 0, 0, 0, 0, 0),
        Vector{NTuple{2,Int}}(undef, 0),
        Vector{NTuple{2,Int}}(undef, 0),
        Matrix{ComplexF64}(undef, 0, 0),
        Vector{NTuple{2,Int}}(undef, 0),
        Vector{NTuple{2,Int}}(undef, 0),
        Matrix{ComplexF64}(undef, 0, 0),
        Matrix{ComplexF64}(undef, 0, 0),
        Vector{Int}(undef, 0),
        Vector{Int}(undef, 0),
        Matrix{ComplexF64}(undef, 0, 0),
        Matrix{ComplexF64}(undef, 0, 0),
        Matrix{ComplexF64}(undef, 0, 0),
        Matrix{ComplexF64}(undef, 0, 0),
        Matrix{ComplexF64}(undef, 0, 0),
        Matrix{ComplexF64}(undef, 0, 0),
        Matrix{ComplexF64}(undef, 0, 0))
end
mutable struct _BSEspinaware_wiTDA_Szn1_buffer <: _BSE_buffer
    ΔE_ud::Vector{Float64} # [nvck_ud]
    WMvvkk_du::Array{ComplexF64,5} # [np_W_dn, nv_up, nv_up, nk, nk]
end
function Base.empty(::Type{_BSEspinaware_wiTDA_Szn1_buffer})
    return _BSEspinaware_wiTDA_Szn1_buffer(
        Vector{Float64}(undef, 0),
        Array{ComplexF64,5}(undef, 0, 0, 0, 0, 0))
end
mutable struct _BSEspinaware_wiTDA_Szp1_buffer <: _BSE_buffer
    ΔE_du::Vector{Float64} # [nvck_du]
    WMvvkk_ud::Array{ComplexF64,5} # [np_W_up, nv_dn, nv_dn, nk, nk]
end
function Base.empty(::Type{_BSEspinaware_wiTDA_Szp1_buffer})
    return _BSEspinaware_wiTDA_Szp1_buffer(
        Vector{Float64}(undef, 0),
        Array{ComplexF64,5}(undef, 0, 0, 0, 0, 0))
end
mutable struct _BSEspinaware_wiTDA_buffer <: _BSE_buffer
    isallocated::Bool
    Kᵈ_ij_up::Vector{NTuple{2,Int}} # [nij_W_up]
    Kᵈ_ij_dn::Vector{NTuple{2,Int}} # [nij_W_dn]
    ekq_pR::Matrix{ComplexF64} #[npR_W, nk]
    Kᵈ_pair_up::Vector{NTuple{2,Int}} # [np_W_up]
    Kᵈ_pair_dn::Vector{NTuple{2,Int}} # [np_W_dn]
    cckk_up::Vector{NTuple{4,Int}} # [nc_up, nc_up, nk, nk]
    cckk_dn::Vector{NTuple{4,Int}} # [nc_dn, nc_dn, nk, nk]
    Mcckk_up::Array{ComplexF64,5} # [np_W_up, nc_up, nc_up, nk, nk]
    Mcckk_dn::Array{ComplexF64,5} # [np_W_dn, nc_dn, nc_dn, nk, nk]
    Sz0::_BSEspinaware_wiTDA_Sz0_buffer
    Szn1::_BSEspinaware_wiTDA_Szn1_buffer
    Szp1::_BSEspinaware_wiTDA_Szp1_buffer
end
function Base.empty(::Type{_BSEspinaware_wiTDA_buffer})
    return _BSEspinaware_wiTDA_buffer(false,
        Vector{NTuple{2,Int}}(undef, 0),
        Vector{NTuple{2,Int}}(undef, 0),
        Matrix{ComplexF64}(undef, 0, 0),
        Vector{NTuple{2,Int}}(undef, 0),
        Vector{NTuple{2,Int}}(undef, 0),
        Vector{NTuple{4,Int}}(undef, 0),
        Vector{NTuple{4,Int}}(undef, 0),
        Array{ComplexF64,5}(undef, 0, 0, 0, 0, 0),
        Array{ComplexF64,5}(undef, 0, 0, 0, 0, 0),
        empty(_BSEspinaware_wiTDA_Sz0_buffer),
        empty(_BSEspinaware_wiTDA_Szn1_buffer),
        empty(_BSEspinaware_wiTDA_Szp1_buffer))
end
function Base.empty!(buffer::_BSEspinaware_wiTDA_buffer)
    buffer.isallocated = false
    buffer.Kᵈ_ij_up = Vector{NTuple{2,Int}}(undef, 0)
    buffer.Kᵈ_ij_dn = Vector{NTuple{2,Int}}(undef, 0)
    buffer.ekq_pR = Matrix{ComplexF64}(undef, 0, 0)
    buffer.Kᵈ_pair_up = Vector{NTuple{2,Int}}(undef, 0)
    buffer.Kᵈ_pair_dn = Vector{NTuple{2,Int}}(undef, 0)
    buffer.cckk_up = Vector{NTuple{4,Int}}(undef, 0)
    buffer.cckk_dn = Vector{NTuple{4,Int}}(undef, 0)
    buffer.Mcckk_up = Array{ComplexF64,5}(undef, 0, 0, 0, 0, 0)
    buffer.Mcckk_dn = Array{ComplexF64,5}(undef, 0, 0, 0, 0, 0)
    buffer.Sz0 = empty(_BSEspinaware_wiTDA_Sz0_buffer)
    buffer.Szn1 = empty(_BSEspinaware_wiTDA_Szn1_buffer)
    buffer.Szp1 = empty(_BSEspinaware_wiTDA_Szp1_buffer)
    return buffer
end
WannierInterActionBase.numorb(bse::BSEspinaware_wiTDA) = length(bse.orb_location_up) + length(bse.orb_location_dn)
function Base.show(io::IO, bse::BSEspinaware_wiTDA)
    print(io, "$(count(bse.period)) dimensinal BSE model with $(numatom(bse)) atoms and $(numorb(bse)) orbitals.")
end
"""
    BSE_wiTDA(kgrid::RedKgrid, TB_up::AbstractTightBindModel, TB_dn::AbstractTightBindModel,
    Kᵈ::PairInterAction, Kˣ::PairInterAction;
    v_up::AbstractVector{<:Integer}=Int[], c_up::AbstractVector{<:Integer}=Int[],
    v_dn::AbstractVector{<:Integer}=Int[], c_dn::AbstractVector{<:Integer}=Int[],
    scissor::Real=0, scissor_up::Real=scissor, scissor_dn::Real=scissor,
    scissor_c_up::Real=scissor_up, scissor_c_dn::Real=scissor_dn, scissor_v_up::Real=0, scissor_v_dn::Real=0)
"""
function BSE_wiTDA(kgrid::RedKgrid, TB_up::AbstractTightBindModel, TB_dn::AbstractTightBindModel,
    Kᵈ::PairInterAction, Kˣ::PairInterAction;
    v_up::AbstractVector{<:Integer}=Int[], c_up::AbstractVector{<:Integer}=Int[],
    v_dn::AbstractVector{<:Integer}=Int[], c_dn::AbstractVector{<:Integer}=Int[],
    scissor::Real=0, scissor_up::Real=scissor, scissor_dn::Real=scissor,
    scissor_c_up::Real=scissor_up, scissor_c_dn::Real=scissor_dn, scissor_v_up::Real=0, scissor_v_dn::Real=0)

    nw_up = numorb(TB_up)
    nw_dn = numorb(TB_dn)
    @assert !iszero(nw_up) "empty up TB is unacceptable!"
    @assert !iszero(nw_dn) "empty down TB is unacceptable!"
    @assert isempty(v_up) || all(i->1 ≤ i ≤ nw_up, v_up) "All `v_up` must lie in [1, $nw_up]!"
    @assert isempty(c_up) || all(i->1 ≤ i ≤ nw_up, c_up) "All `c_up` must lie in [1, $nw_up]!"
    @assert isempty(v_dn) || all(i->1 ≤ i ≤ nw_dn, v_dn) "All `v_dn` must lie in [1, $nw_dn]!"
    @assert isempty(c_dn) || all(i->1 ≤ i ≤ nw_dn, c_dn) "All `c_dn` must lie in [1, $nw_dn]!"
    @assert length(v_up) + length(v_dn) > 0 "Need at least one valence band total (`v_up` + `v_dn`)"
    @assert length(c_up) + length(c_dn) > 0 "Need at least one conduction band total (`c_up` + `c_dn`)"

    nw = nw_up + nw_dn
    @assert (Kᵈ.nw == nw) || (Kᵈ.nw == nw_up == nw_dn) "Only accept Kᵈ.nw=nw_up+nw_dn or Kᵈ.nw=nw_up=nw_dn!"
    @assert (Kˣ.nw == nw) || (Kˣ.nw == nw_up == nw_dn) "Only accept Kˣ.nw=nw_up+nw_dn or Kˣ.nw=nw_up=nw_dn!"

    # Rgrid
    Rgrid, ndegen = gridindex(Val(:WignerSeitz), TB_up.lattice, kgrid.kgrid_size)

    nk = length(kgrid)
    nR = length(Rgrid)

    v_up = collect(Int, v_up)
    c_up = collect(Int, c_up)
    v_dn = collect(Int, v_dn)
    c_dn = collect(Int, c_dn)

    kindex = collect(1:nk)
    vck_uu = TripleList(v_up, c_up, kindex)
    vck_dd = TripleList(v_dn, c_dn, kindex)
    vck_ud = TripleList(v_up, c_dn, kindex)
    vck_du = TripleList(v_dn, c_up, kindex)
    Rindex = collect(1:nR)
    ijR_uu = TripleList(collect(1:nw_up), collect(1:nw_up), Rindex)
    ijR_dd = TripleList(collect(1:nw_dn), collect(1:nw_dn), Rindex)
    ijR_ud = TripleList(collect(1:nw_up), collect(1:nw_dn), Rindex)
    ijR_du = TripleList(collect(1:nw_dn), collect(1:nw_up), Rindex)

    bandk_up = BAND(kgrid, TB_up.H; vector=true)
    bandk_dn = BAND(kgrid, TB_dn.H; vector=true)
    phase_normalize!(:positive_maxabs, bandk_up)
    phase_normalize!(:positive_maxabs, bandk_dn)

    berryconnection_up = _BSE_berryconnection_vck(TB_up.H, kgrid, bandk_up, TB_up.lattice, TB_up.orb_location, vck_uu)
    berryconnection_dn = _BSE_berryconnection_vck(TB_dn.H, kgrid, bandk_dn, TB_up.lattice, TB_dn.orb_location, vck_dd)

    Nvck_Sz0 = length(vck_uu) + length(vck_dd)
    Nvck_Szn1 = length(vck_ud)
    Nvck_Szp1 = length(vck_du)
    Sz = Set{Int}()
    !iszero(Nvck_Sz0) && push!(Sz, 0)
    !iszero(Nvck_Szn1) && push!(Sz, -1)
    !iszero(Nvck_Szp1) && push!(Sz, 1)

    buffer = empty(_BSEspinaware_wiTDA_buffer)
    return BSEspinaware_wiTDA{typeof(TB_up.H),typeof(TB_dn.H),typeof(Kᵈ),typeof(Kˣ),typeof(buffer)}(
        TB_up.lattice, TB_up.atom_name, TB_up.atom_location,
        TB_up.orb_name, TB_up.orb_location, TB_dn.orb_name, TB_dn.orb_location,
        TB_up.period, TB_up.H, TB_dn.H, scissor_c_up, scissor_v_up, scissor_c_dn, scissor_v_dn,
        Kᵈ, Kˣ, kgrid, Rgrid, ndegen, vck_uu, vck_dd, vck_ud, vck_du, ijR_uu, ijR_dd, ijR_ud, ijR_du,
        bandk_up, bandk_dn, Sz, berryconnection_up, berryconnection_dn, buffer)
end
function allocate!(bse::BSEspinaware_wiTDA; force::Bool=false, storage::Symbol=:memory,
    Mmap_Mv::Bool=true, Mmap_Mc::Bool=false, Mmap_Mvck::Bool=false, Mmap_VMvck::Bool=false,
    Mmapfolder::AbstractString="./BSE.Mmap")

    if !force && isallocated(bse)
        return bse
    elseif force && isallocated(bse)
        release!(bse)
    else
        GC.gc()
    end

    Nvck_uu = length(bse.vck_uu)
    Nvck_dd = length(bse.vck_dd)
    Nvck_ud = length(bse.vck_ud)
    Nvck_du = length(bse.vck_du)
    v_up = values(bse.vck_uu)[1]
    c_up = values(bse.vck_uu)[2]
    v_dn = values(bse.vck_dd)[1]
    c_dn = values(bse.vck_dd)[2]
    kgrid = bse.kgrid
    nv_up = length(v_up)
    nc_up = length(c_up)
    nv_dn = length(v_dn)
    nc_dn = length(c_dn)
    nk = length(kgrid)
    Kᵈ = bse.Kᵈ
    Kˣ = bse.Kˣ

    nw_up = numorb(bse.eleH_up)
    nw_dn = numorb(bse.eleH_dn)
    nw = nw_up + nw_dn
    if Kᵈ.nw == nw
        Kᵈ_pair_upindex = findall(p -> (p[1] ≤ nw_up && p[2] ≤ nw_up), Kᵈ.pair)
        Kᵈ_pair_dnindex = findall(p -> (p[1] > nw_up && p[2] > nw_up), Kᵈ.pair)
    else
        Kᵈ_pair_upindex = 1:Kᵈ.np
        Kᵈ_pair_dnindex = 1:Kᵈ.np
    end
    np_W = Kᵈ.np
    np_W_up = length(Kᵈ_pair_upindex)
    np_W_dn = length(Kᵈ_pair_dnindex)
    if Kˣ.nw == nw
        Kˣ_pair_upindex = findall(p -> (p[1] ≤ nw_up && p[2] ≤ nw_up), Kˣ.pair)
        Kˣ_pair_dnindex = findall(p -> (p[1] > nw_up && p[2] > nw_up), Kˣ.pair)
    else
        Kˣ_pair_upindex = 1:Kˣ.np
        Kˣ_pair_dnindex = 1:Kˣ.np
    end
    np_V = Kˣ.np
    np_V_up = length(Kˣ_pair_upindex)
    np_V_dn = length(Kˣ_pair_dnindex)

    # report memory
    if !iszero(Nvck_uu)
        WMvvkk_uu_size = format_memory(sizeof(ComplexF64) * np_W_up * nv_up * nv_up * nk * nk)
        println("WMvvkk_uu: ComplexF64, dims=($np_W_up, $nv_up, $nv_up, $nk, $nk), total=$WMvvkk_uu_size")
    end
    if !iszero(Nvck_dd)
        WMvvkk_dd_size = format_memory(sizeof(ComplexF64) * np_W_dn * nv_dn * nv_dn * nk * nk)
        println("WMvvkk_dd: ComplexF64, dims=($np_W_dn, $nv_dn, $nv_dn, $nk, $nk), total=$WMvvkk_dd_size")
    end
    if Kᵈ.nw == nw
        if !iszero(Nvck_ud)
            WMvvkk_du_size = format_memory(sizeof(ComplexF64) * np_W_dn * nv_up * nv_up * nk * nk)
            println("WMvvkk_du: ComplexF64, dims=($np_W_dn, $nv_up, $nv_up, $nk, $nk), total=$WMvvkk_du_size")
        end
        if !iszero(Nvck_du)
            WMvvkk_ud_size = format_memory(sizeof(ComplexF64) * np_W_up * nv_dn * nv_dn * nk * nk)
            println("WMvvkk_ud: ComplexF64, dims=($np_W_up, $nv_dn, $nv_dn, $nk, $nk), total=$WMvvkk_ud_size")
        end
    end
    if !iszero(nc_up)
        Mcckk_up_size = format_memory(sizeof(ComplexF64) * np_W_up * nc_up * nc_up * nk * nk)
        println("Mcckk_up: ComplexF64, dims=($np_W_up, $nc_up, $nc_up, $nk, $nk), total=$Mcckk_up_size")
    end
    if !iszero(nc_dn)
        Mcckk_dn_size = format_memory(sizeof(ComplexF64) * np_W_dn * nc_dn * nc_dn * nk * nk)
        println("Mcckk_dn: ComplexF64, dims=($np_W_dn, $nc_dn, $nc_dn, $nk, $nk), total=$Mcckk_dn_size")
    end
    if !iszero(Nvck_uu)
        Mvck_uu_size = format_memory(sizeof(ComplexF64) * np_V_up * Nvck_uu)
        println("Mvck_uu: ComplexF64, dims=($np_V_up, $Nvck_uu), total=$Mvck_uu_size")
        VMvck_uu_size = format_memory(sizeof(ComplexF64) * np_V_up * Nvck_uu)
        println("VMvck_uu: ComplexF64, dims=($np_V_up, $Nvck_uu), total=$VMvck_uu_size")
    end
    if !iszero(Nvck_dd)
        Mvck_dd_size = format_memory(sizeof(ComplexF64) * np_V_dn * Nvck_dd)
        println("Mvck_dd: ComplexF64, dims=($np_V_dn, $Nvck_dd), total=$Mvck_dd_size")
        VMvck_dd_size = format_memory(sizeof(ComplexF64) * np_V_dn * Nvck_dd)
        println("VMvck_dd: ComplexF64, dims=($np_V_dn, $Nvck_dd), total=$VMvck_dd_size")
    end
    if Kˣ.nw == nw && !iszero(Nvck_dd)
        VMvck_ud_size = format_memory(sizeof(ComplexF64) * np_V_up * Nvck_dd)
        println("VMvck_ud: ComplexF64, dims=($np_V_up, $Nvck_dd), total=$VMvck_ud_size")
    end
    flush(stdout)

    # allocate
    if storage === :Mmap
        rm(Mmapfolder; recursive=true, force=true)
        mkpath(Mmapfolder)
    end
    if Kᵈ.nw == nw
        if iszero(Nvck_uu)
            WMvvkk_uu = Array{ComplexF64,5}(undef, 0, nv_up, nv_up, nk, nk)
        else
            if storage === :memory
                WMvvkk_uu = Array{ComplexF64,5}(undef, np_W_up, nv_up, nv_up, nk, nk)
            elseif storage === :Mmap
                WMvvkk_uu = Mmap_Mv && !iszero(np_W_up) && !iszero(nv_up) ?
                            Mmap.mmap(joinpath(Mmapfolder, "WMv_uu"), Array{ComplexF64,5}, (np_W_up, nv_up, nv_up, nk, nk)) :
                            Array{ComplexF64,5}(undef, np_W_up, nv_up, nv_up, nk, nk)
            end
        end
        if iszero(Nvck_dd)
            WMvvkk_dd = Array{ComplexF64,5}(undef, 0, nv_dn, nv_dn, nk, nk)
        else
            if storage === :memory
                WMvvkk_dd = Array{ComplexF64,5}(undef, np_W_dn, nv_dn, nv_dn, nk, nk)
            elseif storage === :Mmap
                WMvvkk_dd = Mmap_Mv && !iszero(np_W_dn) && !iszero(nv_dn) ?
                            Mmap.mmap(joinpath(Mmapfolder, "WMv_dd"), Array{ComplexF64,5}, (np_W_dn, nv_dn, nv_dn, nk, nk)) :
                            Array{ComplexF64,5}(undef, np_W_dn, nv_dn, nv_dn, nk, nk)
            end
        end
        if iszero(Nvck_ud)
            WMvvkk_du = Array{ComplexF64,5}(undef, 0, nv_up, nv_up, nk, nk)
        else
            if storage === :memory
                WMvvkk_du = Array{ComplexF64,5}(undef, np_W_dn, nv_up, nv_up, nk, nk)
            elseif storage === :Mmap
                WMvvkk_du = Mmap_Mv && !iszero(np_W_dn) && !iszero(nv_up) ?
                            Mmap.mmap(joinpath(Mmapfolder, "WMv_du"), Array{ComplexF64,5}, (np_W_dn, nv_up, nv_up, nk, nk)) :
                            Array{ComplexF64,5}(undef, np_W_dn, nv_up, nv_up, nk, nk)
            end
        end
        if iszero(Nvck_du)
            WMvvkk_ud = Array{ComplexF64,5}(undef, 0, nv_dn, nv_dn, nk, nk)
        else
            if storage === :memory
                WMvvkk_ud = Array{ComplexF64,5}(undef, np_W_up, nv_dn, nv_dn, nk, nk)
            elseif storage === :Mmap
                WMvvkk_ud = Mmap_Mv && !iszero(np_W_up) && !iszero(nv_dn) ?
                            Mmap.mmap(joinpath(Mmapfolder, "WMv_ud"), Array{ComplexF64,5}, (np_W_up, nv_dn, nv_dn, nk, nk)) :
                            Array{ComplexF64,5}(undef, np_W_up, nv_dn, nv_dn, nk, nk)
            end
        end
    else
        if storage === :memory
            WMvvkk_uu = Array{ComplexF64,5}(undef, np_W_up, nv_up, nv_up, nk, nk)
        elseif storage === :Mmap
            WMvvkk_uu = Mmap_Mv && !iszero(np_W_up) && !iszero(nv_up) ?
                        Mmap.mmap(joinpath(Mmapfolder, "WMv_uu"), Array{ComplexF64,5}, (np_W_up, nv_up, nv_up, nk, nk)) :
                        Array{ComplexF64,5}(undef, np_W_up, nv_up, nv_up, nk, nk)
        end
        if storage === :memory
            WMvvkk_dd = Array{ComplexF64,5}(undef, np_W_dn, nv_dn, nv_dn, nk, nk)
        elseif storage === :Mmap
            WMvvkk_dd = Mmap_Mv && !iszero(np_W_dn) && !iszero(nv_dn) ?
                        Mmap.mmap(joinpath(Mmapfolder, "WMv_dd"), Array{ComplexF64,5}, (np_W_dn, nv_dn, nv_dn, nk, nk)) :
                        Array{ComplexF64,5}(undef, np_W_dn, nv_dn, nv_dn, nk, nk)
        end
        WMvvkk_du = WMvvkk_uu
        WMvvkk_ud = WMvvkk_dd
    end
    if storage === :memory
        Mcckk_up = Array{ComplexF64,5}(undef, np_W_up, nc_up, nc_up, nk, nk)
        Mcckk_dn = Array{ComplexF64,5}(undef, np_W_dn, nc_dn, nc_dn, nk, nk)
        Mvck_uu = Matrix{ComplexF64}(undef, np_V_up, Nvck_uu)
        Mvck_dd = Matrix{ComplexF64}(undef, np_V_dn, Nvck_dd)
        VMvck_uu = Matrix{ComplexF64}(undef, np_V_up, Nvck_uu)
        VMvck_dd = Matrix{ComplexF64}(undef, np_V_dn, Nvck_dd)
        if Kˣ.nw == nw
            VMvck_ud = Matrix{ComplexF64}(undef, np_V_up, Nvck_dd)
        else
            VMvck_ud = VMvck_dd
        end
    elseif storage === :Mmap
        Mcckk_up = Mmap_Mc && !iszero(np_W_up) && !iszero(nc_up) ?
                   Mmap.mmap(joinpath(Mmapfolder, "Mc_up"), Array{ComplexF64,5}, (np_W_up, nc_up, nc_up, nk, nk)) :
                   Array{ComplexF64,5}(undef, np_W_up, nc_up, nc_up, nk, nk)
        Mcckk_dn = Mmap_Mc && !iszero(np_W_dn) && !iszero(nc_dn) ?
                   Mmap.mmap(joinpath(Mmapfolder, "Mc_dn"), Array{ComplexF64,5}, (np_W_dn, nc_dn, nc_dn, nk, nk)) :
                   Array{ComplexF64,5}(undef, np_W_dn, nc_dn, nc_dn, nk, nk)
        Mvck_uu = Mmap_Mvck && !iszero(np_V_up) && !iszero(Nvck_uu) ?
                  Mmap.mmap(joinpath(Mmapfolder, "Mvck_uu"), Matrix{ComplexF64}, (np_V_up, Nvck_uu)) :
                  Matrix{ComplexF64}(undef, np_V_up, Nvck_uu)
        Mvck_dd = Mmap_Mvck && !iszero(np_V_dn) && !iszero(Nvck_dd) ?
                  Mmap.mmap(joinpath(Mmapfolder, "Mvck_dd"), Matrix{ComplexF64}, (np_V_dn, Nvck_dd)) :
                  Matrix{ComplexF64}(undef, np_V_dn, Nvck_dd)
        VMvck_uu = Mmap_VMvck && !iszero(np_V_up) && !iszero(Nvck_uu) ?
                   Mmap.mmap(joinpath(Mmapfolder, "VMvck_uu"), Matrix{ComplexF64}, (np_V_up, Nvck_uu)) :
                   Matrix{ComplexF64}(undef, np_V_up, Nvck_uu)
        VMvck_dd = Mmap_VMvck && !iszero(np_V_dn) && !iszero(Nvck_dd) ?
                   Mmap.mmap(joinpath(Mmapfolder, "VMvck_dd"), Matrix{ComplexF64}, (np_V_dn, Nvck_dd)) :
                   Matrix{ComplexF64}(undef, np_V_dn, Nvck_dd)
        if Kˣ.nw == nw
            VMvck_ud = Mmap_VMvck && !iszero(np_V_up) && !iszero(Nvck_dd) ?
                       Mmap.mmap(joinpath(Mmapfolder, "VMvck_ud"), Matrix{ComplexF64}, (np_V_up, Nvck_dd)) :
                       Matrix{ComplexF64}(undef, np_V_up, Nvck_dd)
        else
            VMvck_ud = VMvck_dd
        end
    end

    if Kᵈ.nw == nw
        Kᵈ_pair_up = Kᵈ.pair[Kᵈ_pair_upindex]
        Kᵈ_pair_dn = map(Kᵈ_pair_dnindex) do ip
            (i, j, ipR) = Kᵈ.pair[ip]
            return (i - nw_up, j - nw_up, ipR)
        end
    else
        Kᵈ_pair_up = Kᵈ.pair
        Kᵈ_pair_dn = Kᵈ.pair
    end
    Kᵈ_ij_up = unique(map(p->(p[1], p[2]), Kᵈ_pair_up))
    Kᵈ_ij_up_idx = Dict{NTuple{2,Int},Int}(ij=>i for (i, ij) in enumerate(Kᵈ_ij_up))
    Kᵈ_pair_up = map(p->(Kᵈ_ij_up_idx[(p[1], p[2])], p[3]), Kᵈ_pair_up)
    Kᵈ_ij_dn = unique(map(p->(p[1], p[2]), Kᵈ_pair_dn))
    Kᵈ_ij_dn_idx = Dict{NTuple{2,Int},Int}(ij=>i for (i, ij) in enumerate(Kᵈ_ij_dn))
    Kᵈ_pair_dn = map(p->(Kᵈ_ij_dn_idx[(p[1], p[2])], p[3]), Kᵈ_pair_dn)

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
    vvkk_up = collect(Iterators.product(1:nv_up, 1:nv_up, 1:nk, 1:nk))
    vvkk_up_chunk = _split_tasks(vvkk_up, nt)
    vvkk_dn = collect(Iterators.product(1:nv_dn, 1:nv_dn, 1:nk, 1:nk))
    vvkk_dn_chunk = _split_tasks(vvkk_dn, nt)
    _BSEspinaware_WMvvkk!(Val(Kᵈ.nw == nw), WMvvkk_uu, WMvvkk_dd, WMvvkk_ud, WMvvkk_du,
        vvkk_up_chunk, vvkk_dn_chunk, bse.bandk_up, bse.bandk_dn, v_up, v_dn,
        ek_pR, Kᵈ_ij_up, Kᵈ_pair_up, Kᵈ_ij_dn, Kᵈ_pair_dn, Wk, kgrid_minusmap,
        Kᵈ_pair_upindex, Kᵈ_pair_dnindex, Nvck_uu, Nvck_dd, Nvck_ud, Nvck_du)

    LinearAlgebra.BLAS.set_num_threads(nt)

    if Kˣ.nw == nw
        Kˣ_pair_up = Kˣ.pair[Kˣ_pair_upindex]
        Kˣ_pair_dn = map(Kˣ_pair_dnindex) do ip
            (i, j, ipR) = Kˣ.pair[ip]
            return (i - nw_up, j - nw_up, ipR)
        end
    else
        Kˣ_pair_up = Kˣ.pair
        Kˣ_pair_dn = Kˣ.pair
    end
    Kˣ_ij_up = unique(map(p->(p[1], p[2]), Kˣ_pair_up))
    Kˣ_ij_up_idx = Dict{NTuple{2,Int},Int}(ij=>i for (i, ij) in enumerate(Kˣ_ij_up))
    Kˣ_pair_up = map(p->(Kˣ_ij_up_idx[(p[1], p[2])], p[3]), Kˣ_pair_up)
    Kˣ_ij_dn = unique(map(p->(p[1], p[2]), Kˣ_pair_dn))
    Kˣ_ij_dn_idx = Dict{NTuple{2,Int},Int}(ij=>i for (i, ij) in enumerate(Kˣ_ij_dn))
    Kˣ_pair_dn = map(p->(Kˣ_ij_dn_idx[(p[1], p[2])], p[3]), Kˣ_pair_dn)

    buffer = bse.buffer
    buffer.isallocated = true
    buffer.Kᵈ_ij_up = Kᵈ_ij_up
    buffer.Kᵈ_ij_dn = Kᵈ_ij_dn
    buffer.ekq_pR = Matrix{ComplexF64}(undef, length(Kᵈ.pR), nk)
    buffer.Kᵈ_pair_up = Kᵈ_pair_up
    buffer.Kᵈ_pair_dn = Kᵈ_pair_dn
    buffer.cckk_up = vec(collect(Iterators.product(1:nc_up, 1:nc_up, 1:nk, 1:nk)))
    buffer.cckk_dn = vec(collect(Iterators.product(1:nc_dn, 1:nc_dn, 1:nk, 1:nk)))
    buffer.Mcckk_up = Mcckk_up
    buffer.Mcckk_dn = Mcckk_dn
    Sz0 = buffer.Sz0
    Szn1 = buffer.Szn1
    Szp1 = buffer.Szp1

    Sz0.ΔE_uu = Vector{Float64}(undef, Nvck_uu)
    Sz0.ΔE_dd = Vector{Float64}(undef, Nvck_dd)
    Sz0.WMvvkk_uu = WMvvkk_uu
    Sz0.WMvvkk_dd = WMvvkk_dd
    Sz0.Kˣ_ij_up = Kˣ_ij_up
    Sz0.Kˣ_ij_dn = Kˣ_ij_dn
    Sz0.ek_pR = [cispi(-2 * (k ⋅ R)) for R in Kˣ.pR, k in kgrid]
    Sz0.Kˣ_pair_up = Kˣ_pair_up
    Sz0.Kˣ_pair_dn = Kˣ_pair_dn
    Sz0.Mvck_uu = Mvck_uu
    Sz0.Mvck_dd = Mvck_dd
    Sz0.Kˣ_pair_upindex = Kˣ_pair_upindex
    Sz0.Kˣ_pair_dnindex = Kˣ_pair_dnindex
    Sz0.Vq = Matrix{ComplexF64}(undef, np_V, np_V)
    Sz0.Vq_uu = Matrix{ComplexF64}(undef, np_V_up, np_V_up)
    Sz0.Vq_dd = Matrix{ComplexF64}(undef, np_V_dn, np_V_dn)
    Sz0.Vq_ud = Matrix{ComplexF64}(undef, np_V_up, np_V_dn)
    Sz0.VMvck_uu = VMvck_uu
    Sz0.VMvck_dd = VMvck_dd
    Sz0.VMvck_ud = VMvck_ud

    Szn1.ΔE_ud = Vector{Float64}(undef, Nvck_ud)
    Szn1.WMvvkk_du = WMvvkk_du

    Szp1.ΔE_du = Vector{Float64}(undef, Nvck_du)
    Szp1.WMvvkk_ud = WMvvkk_ud

    return bse
end
function _BSEspinaware_WMvvkk!(::Val{true}, WMvvkk_uu, WMvvkk_dd, WMvvkk_ud, WMvvkk_du,
    vvkk_up_chunk, vvkk_dn_chunk, bandk_up, bandk_dn, v_up, v_dn,
    ek_pR, Kᵈ_ij_up, Kᵈ_pair_up, Kᵈ_ij_dn, Kᵈ_pair_dn, Wk, kgrid_minusmap,
    Kᵈ_pair_upindex, Kᵈ_pair_dnindex, Nvck_uu, Nvck_dd, Nvck_ud, Nvck_du)
    nv_up = length(v_up)
    nv_dn = length(v_dn)
    nij_W_up = length(Kᵈ_ij_up)
    np_W_up = length(Kᵈ_pair_up)
    nij_W_dn = length(Kᵈ_ij_dn)
    np_W_dn = length(Kᵈ_pair_dn)
    nk = size(Wk, 3)
    Wk_uu = iszero(Nvck_uu) ? Array{ComplexF64}(undef, 0, np_W_up, nk) : Wk[Kᵈ_pair_upindex, Kᵈ_pair_upindex, :]
    Wk_dd = iszero(Nvck_dd) ? Array{ComplexF64}(undef, 0, np_W_dn, nk) : Wk[Kᵈ_pair_dnindex, Kᵈ_pair_dnindex, :]
    Wk_du = iszero(Nvck_ud) ? Array{ComplexF64}(undef, 0, np_W_up, nk) : Wk[Kᵈ_pair_dnindex, Kᵈ_pair_upindex, :]
    Wk_ud = iszero(Nvck_du) ? Array{ComplexF64}(undef, 0, np_W_dn, nk) : Wk[Kᵈ_pair_upindex, Kᵈ_pair_dnindex, :]
    if !iszero(nv_up)
        Threads.@threads for vvkks in vvkk_up_chunk
            buffer_nij_W_up = Vector{ComplexF64}(undef, nij_W_up)
            buffer_np_W_up = Vector{ComplexF64}(undef, np_W_up)
            @inbounds for (iv, iv′, ik, ik′) in vvkks
                ψ = view(bandk_up[ik].vectors, :, v_up[iv])
                ψ′ = view(bandk_up[ik′].vectors, :, v_up[iv′])
                @simd for ij in Base.OneTo(nij_W_up)
                    (i, j) = Kᵈ_ij_up[ij]
                    buffer_nij_W_up[ij] = ψ[i] * conj(ψ′[j])
                end
                ek_pR_k′ = view(ek_pR, :, ik′)
                @simd for ip in Base.OneTo(np_W_up)
                    (ij, ipR) = Kᵈ_pair_up[ip]
                    buffer_np_W_up[ip] = buffer_nij_W_up[ij] * ek_pR_k′[ipR]
                end
                kk′ = kgrid_minusmap[ik, ik′]
                W_kk′_uu = view(Wk_uu, :, :, kk′)
                W_kk′_du = view(Wk_du, :, :, kk′)
                WMvvkk_uu_view = view(WMvvkk_uu, :, iv, iv′, ik, ik′)
                WMvvkk_du_view = view(WMvvkk_du, :, iv, iv′, ik, ik′)
                mul!(WMvvkk_uu_view, W_kk′_uu, buffer_np_W_up)
                mul!(WMvvkk_du_view, W_kk′_du, buffer_np_W_up)
            end
        end
    end
    if !iszero(nv_dn)
        Threads.@threads for vvkks in vvkk_dn_chunk
            buffer_nij_W_dn = Vector{ComplexF64}(undef, nij_W_dn)
            buffer_np_W_dn = Vector{ComplexF64}(undef, np_W_dn)
            @inbounds for (iv, iv′, ik, ik′) in vvkks
                ψ = view(bandk_dn[ik].vectors, :, v_dn[iv])
                ψ′ = view(bandk_dn[ik′].vectors, :, v_dn[iv′])
                @simd for ij in Base.OneTo(nij_W_dn)
                    (i, j) = Kᵈ_ij_dn[ij]
                    buffer_nij_W_dn[ij] = ψ[i] * conj(ψ′[j])
                end
                ek_pR_k′ = view(ek_pR, :, ik′)
                @simd for ip in Base.OneTo(np_W_dn)
                    (ij, ipR) = Kᵈ_pair_dn[ip]
                    buffer_np_W_dn[ip] = buffer_nij_W_dn[ij] * ek_pR_k′[ipR]
                end
                kk′ = kgrid_minusmap[ik, ik′]
                W_kk′_dd = view(Wk_dd, :, :, kk′)
                W_kk′_ud = view(Wk_ud, :, :, kk′)
                WMvvkk_dd_view = view(WMvvkk_dd, :, iv, iv′, ik, ik′)
                WMvvkk_ud_view = view(WMvvkk_ud, :, iv, iv′, ik, ik′)
                mul!(WMvvkk_dd_view, W_kk′_dd, buffer_np_W_dn)
                mul!(WMvvkk_ud_view, W_kk′_ud, buffer_np_W_dn)
            end
        end
    end
end
function _BSEspinaware_WMvvkk!(::Val{false}, WMvvkk_uu, WMvvkk_dd, WMvvkk_ud, WMvvkk_du,
    vvkk_up_chunk, vvkk_dn_chunk, bandk_up, bandk_dn, v_up, v_dn,
    ek_pR, Kᵈ_ij_up, Kᵈ_pair_up, Kᵈ_ij_dn, Kᵈ_pair_dn, Wk, kgrid_minusmap,
    Kᵈ_pair_upindex, Kᵈ_pair_dnindex, Nvck_uu, Nvck_dd, Nvck_ud, Nvck_du)
    nv_up = length(v_up)
    nv_dn = length(v_dn)
    nij_W_up = length(Kᵈ_ij_up)
    np_W_up = length(Kᵈ_pair_up)
    nij_W_dn = length(Kᵈ_ij_dn)
    np_W_dn = length(Kᵈ_pair_dn)
    if !iszero(nv_up) && !iszero(np_W_up)
        Threads.@threads for vvkks in vvkk_up_chunk
            buffer_nij_W_up = Vector{ComplexF64}(undef, nij_W_up)
            buffer_np_W_up = Vector{ComplexF64}(undef, np_W_up)
            @inbounds for (iv, iv′, ik, ik′) in vvkks
                ψ = view(bandk_up[ik].vectors, :, v_up[iv])
                ψ′ = view(bandk_up[ik′].vectors, :, v_up[iv′])
                @simd for ij in Base.OneTo(nij_W_up)
                    (i, j) = Kᵈ_ij_up[ij]
                    buffer_nij_W_up[ij] = ψ[i] * conj(ψ′[j])
                end
                ek_pR_k′ = view(ek_pR, :, ik′)
                @simd for ip in Base.OneTo(np_W_up)
                    (ij, ipR) = Kᵈ_pair_up[ip]
                    buffer_np_W_up[ip] = buffer_nij_W_up[ij] * ek_pR_k′[ipR]
                end
                kk′ = kgrid_minusmap[ik, ik′]
                W_kk′ = view(Wk, :, :, kk′)
                WMvvkk_uu_view = view(WMvvkk_uu, :, iv, iv′, ik, ik′)
                mul!(WMvvkk_uu_view, W_kk′, buffer_np_W_up)
            end
        end
    end
    if !iszero(nv_dn) && !iszero(np_W_dn)
        Threads.@threads for vvkks in vvkk_dn_chunk
            buffer_nij_W_dn = Vector{ComplexF64}(undef, nij_W_dn)
            buffer_np_W_dn = Vector{ComplexF64}(undef, np_W_dn)
            @inbounds for (iv, iv′, ik, ik′) in vvkks
                ψ = view(bandk_dn[ik].vectors, :, v_dn[iv])
                ψ′ = view(bandk_dn[ik′].vectors, :, v_dn[iv′])
                @simd for ij in Base.OneTo(nij_W_dn)
                    (i, j) = Kᵈ_ij_dn[ij]
                    buffer_nij_W_dn[ij] = ψ[i] * conj(ψ′[j])
                end
                ek_pR_k′ = view(ek_pR, :, ik′)
                @simd for ip in Base.OneTo(np_W_dn)
                    (ij, ipR) = Kᵈ_pair_dn[ip]
                    buffer_np_W_dn[ip] = buffer_nij_W_dn[ij] * ek_pR_k′[ipR]
                end
                kk′ = kgrid_minusmap[ik, ik′]
                W_kk′ = view(Wk, :, :, kk′)
                WMvvkk_dd_view = view(WMvvkk_dd, :, iv, iv′, ik, ik′)
                mul!(WMvvkk_dd_view, W_kk′, buffer_np_W_dn)
            end
        end
    end
end
function release!(bse::BSEspinaware_wiTDA)
    isallocated(bse) || return bse
    empty!(bse.buffer)
    GC.gc()
    return bse
end
function _BSE_preprocess_q!(bse::BSEspinaware_wiTDA, q::ReducedCoordinates; Sz=bse.Sz, isΓ::Bool)
    isallocated(bse) || allocate!(bse)
    if isΓ
        # q is direction
        q_value = ReducedCoordinates{Int}(0, 0, 0)
    else
        q_value = q
    end
    has_c_up = false
    has_c_dn = false
    Nvck_uu = length(bse.vck_uu)
    Nvck_dd = length(bse.vck_dd)
    if 0 ∈ Sz
        has_c_up = iszero(Nvck_uu) ? has_c_up : true
        has_c_dn = iszero(Nvck_dd) ? has_c_dn : true
    end
    has_c_dn = -1 ∈ Sz ? true : has_c_dn
    has_c_up = 1 ∈ Sz ? true : has_c_up
    # bandkq
    if isΓ
        bandkq_up = has_c_up ? bse.bandk_up : nothing
        bandkq_dn = has_c_dn ? bse.bandk_dn : nothing
    else
        kq = map(k -> k + q, bse.kgrid)
        if has_c_up
            bandkq_up = BAND(kq, bse.eleH_up; vector=true)
            phase_normalize!(:positive_maxabs, bandkq_up)
        else
            bandkq_up = nothing
        end
        if has_c_dn
            bandkq_dn = BAND(kq, bse.eleH_dn; vector=true)
            phase_normalize!(:positive_maxabs, bandkq_dn)
        else
            bandkq_dn = nothing
        end
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
    cckk_up_chunk = _split_tasks(bse.buffer.cckk_up, nt)
    cckk_dn_chunk = _split_tasks(bse.buffer.cckk_dn, nt)
    _BSEspinaware_Mcckk!(bse.buffer.Mcckk_up, bse.buffer.Mcckk_dn, cckk_up_chunk, cckk_dn_chunk,
        has_c_up, has_c_dn, bandkq_up, bandkq_dn, values(bse.vck_uu)[2], values(bse.vck_dd)[2],
        ekq_pR, bse.buffer.Kᵈ_ij_up, bse.buffer.Kᵈ_pair_up, bse.buffer.Kᵈ_ij_dn, bse.buffer.Kᵈ_pair_dn)

    0 ∈ Sz && _BSEspinaware_preprocess_q!(Val(0), bse, bandkq_up, bandkq_dn, q, isΓ)
    -1 ∈ Sz && _BSEspinaware_preprocess_q!(Val(-1), bse, bandkq_dn)
    1 ∈ Sz && _BSEspinaware_preprocess_q!(Val(1), bse, bandkq_up)

    return bandkq_up, bandkq_dn
end
function _BSEspinaware_Mcckk!(Mcckk_up, Mcckk_dn, cckk_up_chunk, cckk_dn_chunk,
    has_c_up, has_c_dn, bandkq_up, bandkq_dn, c_up, c_dn,
    ekq_pR, Kᵈ_ij_up, Kᵈ_pair_up, Kᵈ_ij_dn, Kᵈ_pair_dn)
    nc_up = length(c_up)
    nc_dn = length(c_dn)
    np_W_up = length(Kᵈ_pair_up)
    np_W_dn = length(Kᵈ_pair_dn)
    if has_c_up && !iszero(nc_up) && !iszero(np_W_up)
        nij_W_up = length(Kᵈ_ij_up)
        Threads.@threads for cckks in cckk_up_chunk
            buffer_nij_W_up = Vector{ComplexF64}(undef, nij_W_up)
            @inbounds for (ic, ic′, ik, ik′) in cckks
                ψ = view(bandkq_up[ik].vectors, :, c_up[ic])
                ψ′ = view(bandkq_up[ik′].vectors, :, c_up[ic′])
                @simd for ij in Base.OneTo(nij_W_up)
                    (i, j) = Kᵈ_ij_up[ij]
                    buffer_nij_W_up[ij] = ψ[i] * conj(ψ′[j])
                end
                Mcckk_up_view = view(Mcckk_up, :, ic, ic′, ik, ik′)
                ekq_pR_k′ = view(ekq_pR, :, ik′)
                @simd for ip in Base.OneTo(np_W_up)
                    (ij, ipR) = Kᵈ_pair_up[ip]
                    Mcckk_up_view[ip] = buffer_nij_W_up[ij] * ekq_pR_k′[ipR]
                end
            end
        end
    end
    if has_c_dn && !iszero(nc_dn) && !iszero(np_W_dn)
        nij_W_dn = length(Kᵈ_ij_dn)
        Threads.@threads for cckks in cckk_dn_chunk
            buffer_nij_W_dn = Vector{ComplexF64}(undef, nij_W_dn)
            @inbounds for (ic, ic′, ik, ik′) in cckks
                ψ = view(bandkq_dn[ik].vectors, :, c_dn[ic])
                ψ′ = view(bandkq_dn[ik′].vectors, :, c_dn[ic′])
                @simd for ij in Base.OneTo(nij_W_dn)
                    (i, j) = Kᵈ_ij_dn[ij]
                    buffer_nij_W_dn[ij] = ψ[i] * conj(ψ′[j])
                end
                Mcckk_dn_view = view(Mcckk_dn, :, ic, ic′, ik, ik′)
                ekq_pR_k′ = view(ekq_pR, :, ik′)
                @simd for ip in Base.OneTo(np_W_dn)
                    (ij, ipR) = Kᵈ_pair_dn[ip]
                    Mcckk_dn_view[ip] = buffer_nij_W_dn[ij] * ekq_pR_k′[ipR]
                end
            end
        end
    end
end
function _BSEspinaware_preprocess_q!(::Val{0}, bse::BSEspinaware_wiTDA, bandkq_up, bandkq_dn, q, isΓ)
    bandk_up = bse.bandk_up
    bandk_dn = bse.bandk_dn
    # ΔE
    vck_uu = bse.vck_uu
    ΔE_uu = bse.buffer.Sz0.ΔE_uu
    scissor = bse.scissor_c_up - bse.scissor_v_up
    for (ivck, (v, c, ik)) in enumerate(vck_uu)
        ΔE_uu[ivck] = bandkq_up[ik].values[c] - bandk_up[ik].values[v] + scissor
    end
    vck_dd = bse.vck_dd
    ΔE_dd = bse.buffer.Sz0.ΔE_dd
    scissor = bse.scissor_c_dn - bse.scissor_v_dn
    for (ivck, (v, c, ik)) in enumerate(vck_dd)
        ΔE_dd[ivck] = bandkq_dn[ik].values[c] - bandk_dn[ik].values[v] + scissor
    end

    nt = Threads.nthreads()

    # Mvck
    Mvck_uu = bse.buffer.Sz0.Mvck_uu
    Mvck_dd = bse.buffer.Sz0.Mvck_dd
    ivck_uu_chunk = _split_tasks(length(vck_uu), nt)
    ivck_dd_chunk = _split_tasks(length(vck_dd), nt)
    _BSEspinaware_Mvck!(Mvck_uu, Mvck_dd, ivck_uu_chunk, ivck_dd_chunk, vck_uu, vck_dd,
        bandk_up, bandk_dn, bandkq_up, bandkq_dn, bse.buffer.Sz0.ek_pR,
        bse.buffer.Sz0.Kˣ_ij_up, bse.buffer.Sz0.Kˣ_pair_up,
        bse.buffer.Sz0.Kˣ_ij_dn, bse.buffer.Sz0.Kˣ_pair_dn)

    # Vq
    nw = numorb(bse)
    Kˣ = bse.Kˣ
    Vq = bse.buffer.Sz0.Vq
    Kˣ(Vq, q; isΓ)
    Vq ./= length(bse.kgrid)
    if Kˣ.nw == nw
        Kˣ_pair_upindex = bse.buffer.Sz0.Kˣ_pair_upindex
        Kˣ_pair_dnindex = bse.buffer.Sz0.Kˣ_pair_dnindex
        Vq_uu = bse.buffer.Sz0.Vq_uu
        Vq_dd = bse.buffer.Sz0.Vq_dd
        Vq_ud = bse.buffer.Sz0.Vq_ud
        VMvck_uu = bse.buffer.Sz0.VMvck_uu
        VMvck_dd = bse.buffer.Sz0.VMvck_dd
        VMvck_ud = bse.buffer.Sz0.VMvck_ud
        Vq_uu .= view(Vq, Kˣ_pair_upindex, Kˣ_pair_upindex)
        Vq_dd .= view(Vq, Kˣ_pair_dnindex, Kˣ_pair_dnindex)
        Vq_ud .= view(Vq, Kˣ_pair_upindex, Kˣ_pair_dnindex)
        mul!(VMvck_uu, Vq_uu, Mvck_uu)
        mul!(VMvck_dd, Vq_dd, Mvck_dd)
        mul!(VMvck_ud, Vq_ud, Mvck_dd)
    else
        VMvck_uu = bse.buffer.Sz0.VMvck_uu
        VMvck_dd = bse.buffer.Sz0.VMvck_dd
        mul!(VMvck_uu, Vq, Mvck_uu)
        mul!(VMvck_dd, Vq, Mvck_dd)
    end

    return nothing
end
function _BSEspinaware_Mvck!(Mvck_uu, Mvck_dd, ivck_uu_chunk, ivck_dd_chunk, vck_uu, vck_dd,
    bandk_up, bandk_dn, bandkq_up, bandkq_dn,
    ek_pR, Kˣ_ij_up, Kˣ_pair_up, Kˣ_ij_dn, Kˣ_pair_dn)
    if !isempty(Mvck_uu)
        nij_V_up = length(Kˣ_ij_up)
        np_V_up = length(Kˣ_pair_up)
        Threads.@threads for ivcks in ivck_uu_chunk
            buffer_nij_V_up = Vector{ComplexF64}(undef, nij_V_up)
            @inbounds for ivck in ivcks
                (v, c, ik) = vck_uu[ivck]
                ψc = view(bandkq_up[ik].vectors, :, c)
                ψv = view(bandk_up[ik].vectors, :, v)
                @simd for ij in Base.OneTo(nij_V_up)
                    (i, j) = Kˣ_ij_up[ij]
                    buffer_nij_V_up[ij] = ψc[i] * conj(ψv[j])
                end
                Mvck_uu_view = view(Mvck_uu, :, ivck)
                ek_pR_k = view(ek_pR, :, ik)
                @simd for ip in Base.OneTo(np_V_up)
                    (ij, ipR) = Kˣ_pair_up[ip]
                    Mvck_uu_view[ip] = buffer_nij_V_up[ij] * ek_pR_k[ipR]
                end
            end
        end
    end
    if !isempty(Mvck_dd)
        nij_V_dn = length(Kˣ_ij_dn)
        np_V_dn = length(Kˣ_pair_dn)
        Threads.@threads for ivcks in ivck_dd_chunk
            buffer_nij_V_dn = Vector{ComplexF64}(undef, nij_V_dn)
            @inbounds for ivck in ivcks
                (v, c, ik) = vck_dd[ivck]
                ψc = view(bandkq_dn[ik].vectors, :, c)
                ψv = view(bandk_dn[ik].vectors, :, v)
                @simd for ij in Base.OneTo(nij_V_dn)
                    (i, j) = Kˣ_ij_dn[ij]
                    buffer_nij_V_dn[ij] = ψc[i] * conj(ψv[j])
                end
                Mvck_dd_view = view(Mvck_dd, :, ivck)
                ek_pR_k = view(ek_pR, :, ik)
                @simd for ip in Base.OneTo(np_V_dn)
                    (ij, ipR) = Kˣ_pair_dn[ip]
                    Mvck_dd_view[ip] = buffer_nij_V_dn[ij] * ek_pR_k[ipR]
                end
            end
        end
    end
end
function _BSEspinaware_preprocess_q!(::Val{-1}, bse::BSEspinaware_wiTDA, bandkq_dn)
    bandk_up = bse.bandk_up
    vck_ud = bse.vck_ud
    ΔE_ud = bse.buffer.Szn1.ΔE_ud
    scissor = bse.scissor_c_dn - bse.scissor_v_up
    for (ivck, (v, c, ik)) in enumerate(vck_ud)
        ΔE_ud[ivck] = bandkq_dn[ik].values[c] - bandk_up[ik].values[v] + scissor
    end
    return nothing
end
function _BSEspinaware_preprocess_q!(::Val{1}, bse::BSEspinaware_wiTDA, bandkq_up)
    bandk_dn = bse.bandk_dn
    vck_du = bse.vck_du
    ΔE_du = bse.buffer.Szp1.ΔE_du
    scissor = bse.scissor_c_up - bse.scissor_v_dn
    for (ivck, (v, c, ik)) in enumerate(vck_du)
        ΔE_du[ivck] = bandkq_up[ik].values[c] - bandk_dn[ik].values[v] + scissor
    end
    return nothing
end
function (bse::BSEspinaware_wiTDA)(::Val{:dimension})
    return bse(Val(:dimension_Sz0)), bse(Val(:dimension_Szn1)), bse(Val(:dimension_Szp1))
end
(bse::BSEspinaware_wiTDA)(::Val{:dimension_Sz0}) = length(bse.vck_uu) + length(bse.vck_dd)
(bse::BSEspinaware_wiTDA)(::Val{:dimension_Szn1}) = length(bse.vck_ud)
(bse::BSEspinaware_wiTDA)(::Val{:dimension_Szp1}) = length(bse.vck_du)
# isΓ表示是否精确处理Γ点，这时输入的q表示方向，注意是分数坐标。
function (bse::BSEspinaware_wiTDA)(q::ReducedCoordinates; isΓ::Bool=false)
    (Nvck_Sz0, Nvck_Szn1, Nvck_Szp1) = bse(Val(:dimension))
    (q, isΓ) = _BSE_preprocess_Γ(q, isΓ, bse.period)
    bandkq_up, bandkq_dn = _BSE_preprocess_q!(bse, q; isΓ)
    HSz0 = iszero(Nvck_Sz0) ? nothing :
           _BSE_Hamiltonian!(Val(0), Matrix{ComplexF64}(undef, Nvck_Sz0, Nvck_Sz0), bse, q; isΓ)
    HSzn1 = iszero(Nvck_Szn1) ? nothing :
            _BSE_Hamiltonian!(Val(-1), Matrix{ComplexF64}(undef, Nvck_Szn1, Nvck_Szn1), bse)
    HSzp1 = iszero(Nvck_Szp1) ? nothing :
            _BSE_Hamiltonian!(Val(1), Matrix{ComplexF64}(undef, Nvck_Szp1, Nvck_Szp1), bse)
    return (bandkq_up=bandkq_up, bandkq_dn=bandkq_dn, HSz0=HSz0, HSzn1=HSzn1, HSzp1=HSzp1)
end
function (bse::BSEspinaware_wiTDA)(Sz::Union{Integer,AbstractVector{<:Integer}},
    q::ReducedCoordinates; isΓ::Bool=false)
    if Sz isa Integer
        Sz = [Sz]
    end
    @assert issubset(Sz, (0, -1, 1)) || "Sz only can contain 0, -1 and 1"
    Sz = intersect(bse.Sz, Sz)
    (Nvck_Sz0, Nvck_Szn1, Nvck_Szp1) = bse(Val(:dimension))
    has_Sz0 = 0 ∈ Sz
    has_Szn1 = -1 ∈ Sz
    has_Szp1 = 1 ∈ Sz
    if !has_Sz0 && !has_Szn1 && !has_Szp1
        return (bandkq_up=nothing, bandkq_dn=nothing, HSz0=nothing, HSzn1=nothing, HSzp1=nothing)
    end
    (q, isΓ) = _BSE_preprocess_Γ(q, isΓ, bse.period)
    bandkq_up, bandkq_dn = _BSE_preprocess_q!(bse, q; Sz, isΓ)
    HSz0 = has_Sz0 ? _BSE_Hamiltonian!(Val(0), Matrix{ComplexF64}(undef, Nvck_Sz0, Nvck_Sz0), bse, q; isΓ) : nothing
    HSzn1 = has_Szn1 ? _BSE_Hamiltonian!(Val(-1), Matrix{ComplexF64}(undef, Nvck_Szn1, Nvck_Szn1), bse) : nothing
    HSzp1 = has_Szp1 ? _BSE_Hamiltonian!(Val(1), Matrix{ComplexF64}(undef, Nvck_Szp1, Nvck_Szp1), bse) : nothing
    return (bandkq_up=bandkq_up, bandkq_dn=bandkq_dn, HSz0=HSz0, HSzn1=HSzn1, HSzp1=HSzp1)
end
function (bse::BSEspinaware_wiTDA)(Hlist::NamedTuple, q::ReducedCoordinates; isΓ::Bool=false)
    @assert issubset(keys(Hlist), (:HSz0, :HSzn1, :HSzp1)) || "Hlist only can contain :HSz0, :HSzn1 and :HSzp1"
    (Nvck_Sz0, Nvck_Szn1, Nvck_Szp1) = bse(Val(:dimension))
    HSz0 = get(Hlist, :HSz0, nothing)
    has_Sz0 = !iszero(Nvck_Sz0) && HSz0 isa AbstractMatrix{<:Number}
    if has_Sz0
        size(HSz0) == (Nvck_Sz0, Nvck_Sz0) ||
            throw(DimensionMismatch("HSz0 must be a $Nvck_Sz0×$Nvck_Sz0 matrix, got size $(size(HSz0))"))
    end
    HSzn1 = get(Hlist, :HSzn1, nothing)
    has_Szn1 = !iszero(Nvck_Szn1) && HSzn1 isa AbstractMatrix{<:Number}
    if has_Szn1
        size(HSzn1) == (Nvck_Szn1, Nvck_Szn1) ||
            throw(DimensionMismatch("HSzn1 must be a $Nvck_Szn1×$Nvck_Szn1 matrix, got size $(size(HSzn1))"))
    end
    HSzp1 = get(Hlist, :HSzp1, nothing)
    has_Szp1 = !iszero(Nvck_Szp1) && HSzp1 isa AbstractMatrix{<:Number}
    if has_Szp1
        size(HSzp1) == (Nvck_Szp1, Nvck_Szp1) ||
            throw(DimensionMismatch("HSzp1 must be a $Nvck_Szp1×$Nvck_Szp1 matrix, got size $(size(HSzp1))"))
    end
    if !has_Sz0 && !has_Szn1 && !has_Szp1
        return (bandkq_up=nothing, bandkq_dn=nothing, HSz0=nothing, HSzn1=nothing, HSzp1=nothing)
    end
    Sz = Int[]
    has_Sz0 && push!(Sz, 0)
    has_Szn1 && push!(Sz, -1)
    has_Szp1 && push!(Sz, 1)
    (q, isΓ) = _BSE_preprocess_Γ(q, isΓ, bse.period)
    bandkq_up, bandkq_dn = _BSE_preprocess_q!(bse, q; Sz, isΓ)
    HSz0 = has_Sz0 ? _BSE_Hamiltonian!(Val(0), HSz0, bse, q; isΓ) : nothing
    HSzn1 = has_Szn1 ? _BSE_Hamiltonian!(Val(-1), HSzn1, bse) : nothing
    HSzp1 = has_Szp1 ? _BSE_Hamiltonian!(Val(1), HSzp1, bse) : nothing
    return (bandkq_up=bandkq_up, bandkq_dn=bandkq_dn, HSz0=HSz0, HSzn1=HSzn1, HSzp1=HSzp1)
end
function _BSE_Hamiltonian!(::Val{0}, H0::AbstractMatrix, bse::BSEspinaware_wiTDA, q; isΓ)
    vck_uu = bse.vck_uu
    vck_dd = bse.vck_dd
    Nvck_uu = length(vck_uu)
    Nvck_dd = length(vck_dd)
    Mcckk_up = bse.buffer.Mcckk_up
    Mcckk_dn = bse.buffer.Mcckk_dn
    ΔE_uu = bse.buffer.Sz0.ΔE_uu
    ΔE_dd = bse.buffer.Sz0.ΔE_dd
    WMvvkk_uu = bse.buffer.Sz0.WMvvkk_uu
    WMvvkk_dd = bse.buffer.Sz0.WMvvkk_dd
    Mvck_uu = bse.buffer.Sz0.Mvck_uu
    Mvck_dd = bse.buffer.Sz0.Mvck_dd
    VMvck_uu = bse.buffer.Sz0.VMvck_uu
    VMvck_dd = bse.buffer.Sz0.VMvck_dd
    VMvck_ud = bse.buffer.Sz0.VMvck_ud
    LinearAlgebra.BLAS.set_num_threads(1)
    Threads.@threads for (ivck, ivck′) in UpperTriangularIndices(Nvck_uu)
        (iv, ic, ik) = triple_index(vck_uu, ivck)
        (iv′, ic′, ik′) = triple_index(vck_uu, ivck′)
        Mc = view(Mcckk_up, :, ic, ic′, ik, ik′)
        WMv = view(WMvvkk_uu, :, iv, iv′, ik, ik′)
        Kᵈ = - Mc ⋅ WMv
        M = view(Mvck_uu, :, ivck)
        VM′ = view(VMvck_uu, :, ivck′)
        Kˣ = M ⋅ VM′
        if ivck == ivck′
            H0[ivck, ivck] = ΔE_uu[ivck] + Kᵈ + Kˣ
        else
            H0[ivck, ivck′] = Kᵈ + Kˣ
        end
    end
    Threads.@threads for (ivck, ivck′) in UpperTriangularIndices(Nvck_dd)
        (iv, ic, ik) = triple_index(vck_dd, ivck)
        (iv′, ic′, ik′) = triple_index(vck_dd, ivck′)
        Mc = view(Mcckk_dn, :, ic, ic′, ik, ik′)
        WMv = view(WMvvkk_dd, :, iv, iv′, ik, ik′)
        Kᵈ = - Mc ⋅ WMv
        M = view(Mvck_dd, :, ivck)
        VM′ = view(VMvck_dd, :, ivck′)
        Kˣ = M ⋅ VM′
        if ivck == ivck′
            i = ivck + Nvck_uu
            H0[i, i] = ΔE_dd[ivck] + Kᵈ + Kˣ
        else
            H0[ivck+Nvck_uu, ivck′+Nvck_uu] = Kᵈ + Kˣ
        end
    end
    Threads.@threads for I in CartesianIndices((Nvck_uu, Nvck_dd))
        (ivck, ivck′) = Tuple(I)
        M = view(Mvck_uu, :, ivck)
        VM′ = view(VMvck_ud, :, ivck′)
        Kˣ = M ⋅ VM′
        H0[ivck, ivck′+Nvck_uu] = Kˣ
    end
    # infinite term in Kˣ
    if isΓ && count(bse.period) == 3
        _BSE_KˣΓ!(H0, bse, q)
    end
    LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
    return Hermitian(H0, :U)
end
function _BSE_KˣΓ!(H0, bse::BSEspinaware_wiTDA, q)
    Ω = abs(det(parent(bse.lattice)))
    CoulombScaledivnk = CoulombScale * 4π / Ω / length(bse.kgrid)
    rlattice = reciprocal(bse.lattice)
    qcar = rlattice * q
    qnorm = norm(qcar)
    q̂x = qcar[1] / qnorm
    q̂y = qcar[2] / qnorm
    q̂z = qcar[3] / qnorm
    A_up = bse.berryconnection_up
    A_dn = bse.berryconnection_dn
    Nvck_uu = length(bse.vck_uu)
    Nvck_dd = length(bse.vck_dd)
    Threads.@threads for (ivck, ivck′) in UpperTriangularIndices(Nvck_uu)
        H0[ivck, ivck′] += conj(A_up[1, ivck] * q̂x + A_up[2, ivck] * q̂y + A_up[3, ivck] * q̂z) *
            (A_up[1, ivck′] * q̂x + A_up[2, ivck′] * q̂y + A_up[3, ivck′] * q̂z) * CoulombScaledivnk
    end
    Threads.@threads for (ivck, ivck′) in UpperTriangularIndices(Nvck_dd)
        H0[ivck+Nvck_uu, ivck′+Nvck_uu] += conj(A_dn[1, ivck] * q̂x + A_dn[2, ivck] * q̂y + A_dn[3, ivck] * q̂z) *
            (A_dn[1, ivck′] * q̂x + A_dn[2, ivck′] * q̂y + A_dn[3, ivck′] * q̂z) * CoulombScaledivnk
    end
    Threads.@threads for I in CartesianIndices((Nvck_uu, Nvck_dd))
        (ivck, ivck′) = Tuple(I)
        H0[ivck, ivck′+Nvck_uu] += conj(A_up[1, ivck] * q̂x + A_up[2, ivck] * q̂y + A_up[3, ivck] * q̂z) *
            (A_dn[1, ivck′] * q̂x + A_dn[2, ivck′] * q̂y + A_dn[3, ivck′] * q̂z) * CoulombScaledivnk
    end
    return H0
end
function _BSE_Hamiltonian!(::Val{-1}, Hn1::AbstractMatrix, bse::BSEspinaware_wiTDA)
    vck_ud = bse.vck_ud
    Mcckk_dn = bse.buffer.Mcckk_dn
    ΔE_ud = bse.buffer.Szn1.ΔE_ud
    WMvvkk_du = bse.buffer.Szn1.WMvvkk_du
    LinearAlgebra.BLAS.set_num_threads(1)
    Threads.@threads for (ivck, ivck′) in UpperTriangularIndices(length(vck_ud))
        (iv, ic, ik) = triple_index(vck_ud, ivck)
        (iv′, ic′, ik′) = triple_index(vck_ud, ivck′)
        Mc = view(Mcckk_dn, :, ic, ic′, ik, ik′)
        WMv = view(WMvvkk_du, :, iv, iv′, ik, ik′)
        Kᵈ = - Mc ⋅ WMv
        if ivck == ivck′
            Hn1[ivck, ivck] = ΔE_ud[ivck] + Kᵈ
        else
            Hn1[ivck, ivck′] = Kᵈ
        end
    end
    LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
    return Hermitian(Hn1, :U)
end
function _BSE_Hamiltonian!(::Val{1}, Hp1::AbstractMatrix, bse::BSEspinaware_wiTDA)
    vck_du = bse.vck_du
    Mcckk_up = bse.buffer.Mcckk_up
    ΔE_du = bse.buffer.Szp1.ΔE_du
    WMvvkk_ud = bse.buffer.Szp1.WMvvkk_ud
    LinearAlgebra.BLAS.set_num_threads(1)
    Threads.@threads for (ivck, ivck′) in UpperTriangularIndices(length(vck_du))
        (iv, ic, ik) = triple_index(vck_du, ivck)
        (iv′, ic′, ik′) = triple_index(vck_du, ivck′)
        Mc = view(Mcckk_up, :, ic, ic′, ik, ik′)
        WMv = view(WMvvkk_ud, :, iv, iv′, ik, ik′)
        Kᵈ = - Mc ⋅ WMv
        if ivck == ivck′
            Hp1[ivck, ivck] = ΔE_du[ivck] + Kᵈ
        else
            Hp1[ivck, ivck′] = Kᵈ
        end
    end
    LinearAlgebra.BLAS.set_num_threads(Threads.nthreads())
    return Hermitian(Hp1, :U)
end
