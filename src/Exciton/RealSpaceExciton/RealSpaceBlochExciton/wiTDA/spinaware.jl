"""
    realspace(bse::BSEspinaware_wiTDA, q::AbstractVector{<:Real}, A::AbstractVector{<:Number}, Sz::Integer) -> RealSpaceExciton

Construct the real-space representation of an exciton in the electron-hole Wannier basis.
This method is intended for the spin-polarized cases, where Wannier functions of both spin
species are included.

The exciton wavefunction with momentum `q` (in fractional coordinates, of length 3) and
total spin `Sz` is given by the coefficients `A` in the electron-hole band basis and is
expanded in the product basis of electron and hole Wannier functions. `Sz` takes one of the
values `-1`, `0`, and `1`.

The returned object `ψ` is a subtype of [`RealSpaceExciton`](@ref) and is callable as

    ψ(ih, Rₕ, ie, Rₑ)

where `ih` and `ie` are the integer Wannier-function indices of the hole and the electron,
respectively, and `Rₕ` and `Rₑ` are three-integer vectors specifying the lattice cells
where the corresponding Wannier functions are centered. Since two spin species are present,
both `ih` and `ie` range over `1:(nw_up + nw_dn)`: the indices `1:nw_up` refer to the
spin-up Wannier functions, and `nw_up+1:nw_up+nw_dn` to the spin-down ones.
"""
function WannierInterActionBase.realspace(bse::BSEspinaware_wiTDA, q::AbstractVector{<:Real}, A::AbstractVector{<:Number}, Sz::Integer)
    return realspace(bse, ReducedCoordinates(q), A, Val(Sz))
end
struct RealSpaceBlochExciton_Sz0_wiTDA <: RealSpaceExciton
    q::ReducedCoordinates{Float64}
    nk::Int
    nw_up::Int
    nw_dn::Int
    inv_sqrt_nk::Float64
    kgrid::RedKgrid
    AU_uu::Array{ComplexF64} #∑_{vc} AUU* 
    AU_dd::Array{ComplexF64} #∑_{vc} AUU* 
end
function (ψ::RealSpaceBlochExciton_Sz0_wiTDA)(ih, Rₕ, ie, Rₑ)
    if ih ≤ ψ.nw_up && ie ≤ ψ.nw_up
        phase = cispi(2 * ψ.q ⋅ Rₑ)
        return ψ.inv_sqrt_nk * phase * sum(Base.OneTo(ψ.nk)) do ik
            ψ.AU_uu[ik, ih, ie] * cispi(2 * ψ.kgrid[ik] ⋅ (Rₑ - Rₕ))
        end
    elseif ih > ψ.nw_up && ie > ψ.nw_up
        ih -= ψ.nw_up
        ie -= ψ.nw_up
        phase = cispi(2 * ψ.q ⋅ Rₑ)
        return ψ.inv_sqrt_nk * phase * sum(Base.OneTo(ψ.nk)) do ik
            ψ.AU_dd[ik, ih, ie] * cispi(2 * ψ.kgrid[ik] ⋅ (Rₑ - Rₕ))
        end
    else
        return 0
    end
end
function WannierInterActionBase.realspace(bse::BSEspinaware_wiTDA, q::ReducedCoordinates, A::AbstractVector{<:Number}, ::Val{0})
    vck_uu = bse.vck_uu
    Nvck_uu = length(vck_uu)
    vck_dd = bse.vck_dd
    Nvck_dd = length(vck_dd)
    length(A) == Nvck_uu + Nvck_dd || throw(DimensionMismatch(
        "length of exciton wavefunction should be Nvck=$(Nvck_uu + Nvck_dd), got $(length(A))"))

    kgrid = bse.kgrid
    nk = length(kgrid)
    bandk_up = bse.bandk_up
    bandk_dn = bse.bandk_dn
    kq = map(k -> k + q, bse.kgrid)
    norb_up = numorb(bse.eleH_up)
    norb_dn = numorb(bse.eleH_dn)

    if iszero(Nvck_uu)
        AU_uu = zeros(ComplexF64, nk, norb_up, norb_up) # ik, ih, ie
    else
        bandkq_up = BAND(kq, bse.eleH_up; vector=true)
        phase_normalize!(:positive_maxabs, bandkq_up)
        (nv, nc, _) = vsize(vck_uu)
        Avck_uu = reshape(view(A, 1:Nvck_uu), nv, nc, nk)

        v_up = values(vck_uu)[1]
        c_up = values(vck_uu)[2]
        vc_idx = [(v, c) for v in 1:nv, c in 1:nc]

        AU_uu = Array{ComplexF64}(undef, nk, norb_up, norb_up) # ik, ih, ie
        Threads.@threads for I in CartesianIndices(AU_uu)
            (ik, ih, ie) = Tuple(I)
            Uv = conj.(bandk_up[ik].vectors[ih, v_up])
            Uc = bandkq_up[ik].vectors[ie, c_up]
            AU_uu[I] = sum(vc_idx) do (v, c)
                Avck_uu[v, c, ik] * Uv[v] * Uc[c]
            end
        end
    end

    if iszero(Nvck_dd)
        AU_dd = zeros(ComplexF64, nk, norb_dn, norb_dn) # ik, ih, ie
    else
        bandkq_dn = BAND(kq, bse.eleH_dn; vector=true)
        phase_normalize!(:positive_maxabs, bandkq_dn)
        (nv, nc, _) = vsize(vck_dd)
        Avck_dd = reshape(view(A, (Nvck_uu+1):(Nvck_uu+Nvck_dd)), nv, nc, nk)

        v_dn = values(vck_dd)[1]
        c_dn = values(vck_dd)[2]
        vc_idx = [(v, c) for v in 1:nv, c in 1:nc]

        AU_dd = Array{ComplexF64}(undef, nk, norb_dn, norb_dn) # ik, ih, ie
        Threads.@threads for I in CartesianIndices(AU_dd)
            (ik, ih, ie) = Tuple(I)
            Uv = conj.(bandk_dn[ik].vectors[ih, v_dn])
            Uc = bandkq_dn[ik].vectors[ie, c_dn]
            AU_dd[I] = sum(vc_idx) do (v, c)
                Avck_dd[v, c, ik] * Uv[v] * Uc[c]
            end
        end
    end

    return RealSpaceBlochExciton_Sz0_wiTDA(q, nk, norb_up, norb_dn, 1/sqrt(nk), kgrid, AU_uu, AU_dd)
end
struct RealSpaceBlochExciton_Szp1_wiTDA <: RealSpaceExciton
    q::ReducedCoordinates{Float64}
    nk::Int
    nw_up::Int
    nw_dn::Int
    inv_sqrt_nk::Float64
    kgrid::RedKgrid
    AU_du::Array{ComplexF64} #∑_{vc} AUU* 
end
function (ψ::RealSpaceBlochExciton_Szp1_wiTDA)(ih, Rₕ, ie, Rₑ)
    if ih > ψ.nw_up && ie ≤ ψ.nw_up
        ih -= ψ.nw_up
        phase = cispi(2 * ψ.q ⋅ Rₑ)
        return ψ.inv_sqrt_nk * phase * sum(Base.OneTo(ψ.nk)) do ik
            ψ.AU_du[ik, ih, ie] * cispi(2 * ψ.kgrid[ik] ⋅ (Rₑ - Rₕ))
        end
    else
        return 0
    end
end
function WannierInterActionBase.realspace(bse::BSEspinaware_wiTDA, q::ReducedCoordinates, A::AbstractVector{<:Number}, ::Val{1})
    vck_du = bse.vck_du
    length(A) == length(vck_du) || throw(DimensionMismatch(
        "length of exciton wavefunction should be Nvck=$(length(vck_du)), got $(length(A))"))

    kgrid = bse.kgrid
    bandk_dn = bse.bandk_dn
    kq = map(k -> k + q, bse.kgrid)
    bandkq_up = BAND(kq, bse.eleH_up; vector=true)
    phase_normalize!(:positive_maxabs, bandkq_up)

    (nv, nc, nk) = vsize(vck_du)
    Avck_du = reshape(A, nv, nc, nk)

    v_dn = values(vck_du)[1]
    c_up = values(vck_du)[2]
    vc_idx = [(v, c) for v in 1:nv, c in 1:nc]
    norb_up = numorb(bse.eleH_up)
    norb_dn = numorb(bse.eleH_dn)
    AU_du = Array{ComplexF64}(undef, nk, norb_dn, norb_up) # ik, ih, ie
    Threads.@threads for I in CartesianIndices(AU_du)
        (ik, ih, ie) = Tuple(I)
        Uv = conj.(bandk_dn[ik].vectors[ih, v_dn])
        Uc = bandkq_up[ik].vectors[ie, c_up]
        AU_du[I] = sum(vc_idx) do (v, c)
            Avck_du[v, c, ik] * Uv[v] * Uc[c]
        end
    end

    return RealSpaceBlochExciton_Szp1_wiTDA(q, nk, norb_up, norb_dn, 1/sqrt(nk), kgrid, AU_du)
end
struct RealSpaceBlochExciton_Szn1_wiTDA <: RealSpaceExciton
    q::ReducedCoordinates{Float64}
    nk::Int
    nw_up::Int
    nw_dn::Int
    inv_sqrt_nk::Float64
    kgrid::RedKgrid
    AU_ud::Array{ComplexF64} #∑_{vc} AUU* 
end
function (ψ::RealSpaceBlochExciton_Szn1_wiTDA)(ih, Rₕ, ie, Rₑ)
    if ih ≤ ψ.nw_up && ie > ψ.nw_up
        ie -= ψ.nw_up
        phase = cispi(2 * ψ.q ⋅ Rₑ)
        return ψ.inv_sqrt_nk * phase * sum(Base.OneTo(ψ.nk)) do ik
            ψ.AU_ud[ik, ih, ie] * cispi(2 * ψ.kgrid[ik] ⋅ (Rₑ - Rₕ))
        end
    else
        return 0
    end
end
function WannierInterActionBase.realspace(bse::BSEspinaware_wiTDA, q::ReducedCoordinates, A::AbstractVector{<:Number}, ::Val{-1})
    vck_ud = bse.vck_ud
    length(A) == length(vck_ud) || throw(DimensionMismatch(
        "length of exciton wavefunction should be Nvck=$(length(vck_ud)), got $(length(A))"))

    kgrid = bse.kgrid
    bandk_up = bse.bandk_up
    kq = map(k -> k + q, bse.kgrid)
    bandkq_dn = BAND(kq, bse.eleH_dn; vector=true)
    phase_normalize!(:positive_maxabs, bandkq_dn)

    (nv, nc, nk) = vsize(vck_ud)
    Avck_ud = reshape(A, nv, nc, nk)

    v_up = values(vck_ud)[1]
    c_dn = values(vck_ud)[2]
    vc_idx = [(v, c) for v in 1:nv, c in 1:nc]
    norb_up = numorb(bse.eleH_up)
    norb_dn = numorb(bse.eleH_dn)
    AU_ud = Array{ComplexF64}(undef, nk, norb_up, norb_dn) # ik, ih, ie
    Threads.@threads for I in CartesianIndices(AU_ud)
        (ik, ih, ie) = Tuple(I)
        Uv = conj.(bandk_up[ik].vectors[ih, v_up])
        Uc = bandkq_dn[ik].vectors[ie, c_dn]
        AU_ud[I] = sum(vc_idx) do (v, c)
            Avck_ud[v, c, ik] * Uv[v] * Uc[c]
        end
    end

    return RealSpaceBlochExciton_Szn1_wiTDA(q, nk, norb_up, norb_dn, 1/sqrt(nk), kgrid, AU_ud)
end
