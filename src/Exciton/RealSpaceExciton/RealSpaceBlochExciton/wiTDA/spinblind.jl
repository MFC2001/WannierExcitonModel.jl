struct RealSpaceBlochExciton_spinblind_wiTDA <: RealSpaceExciton
    q::ReducedCoordinates{Float64}
    nk::Int
    inv_sqrt_nk::Float64
    kgrid::RedKgrid
    AU::Array{ComplexF64} #∑_{vc} AUU* 
end
function (ψ::RealSpaceBlochExciton_spinblind_wiTDA)(ih, Rₕ, ie, Rₑ)
    inv_sqrt_nk_phase = ψ.inv_sqrt_nk * cispi(2 * ψ.q ⋅ Rₑ)
    ΔR = Rₑ - Rₕ
    return inv_sqrt_nk_phase * sum(Base.OneTo(ψ.nk)) do ik
        ψ.AU[ik, ih, ie] * cispi(2 * ψ.kgrid[ik] ⋅ ΔR)
    end
end
"""
    realspace(bse::BSEspinblind_wiTDA, q::AbstractVector{<:Real}, A::AbstractVector{<:Number}) -> RealSpaceExciton

Construct the real-space representation of an exciton in the electron-hole Wannier basis.
This method is intended for the spin-degenerate and noncollinear cases.

The exciton wavefunction with momentum `q` (in fractional coordinates, of length 3)
given by the coefficients `A` in the electron-hole band basis is
expanded in the product basis of electron and hole Wannier functions. 
The returned object `ψ` is a subtype of [`RealSpaceExciton`](@ref) and is callable as

    ψ(ih, Rₕ, ie, Rₑ)

where `ih` and `ie` are the integer Wannier-function indices of the hole and the electron,
respectively, and `Rₕ` and `Rₑ` are three-integer vectors specifying the lattice cells
where the corresponding Wannier functions are centered.
"""
function WannierInterActionBase.realspace(bse::BSEspinblind_wiTDA, q::AbstractVector{<:Real}, A::AbstractVector{<:Number})
    vck = bse.vck
    length(A) == length(vck) || throw(DimensionMismatch(
        "length of exciton wavefunction should be Nvck=$(length(vck)), got $(length(A))"))

    kgrid = bse.kgrid
    bandk = bse.bandk
    q = ReducedCoordinates(q)
    bandkq = BAND(map(k -> k + q, kgrid), bse.eleH; vector=true)
    phase_normalize!(:positive_maxabs, bandkq)

    v = values(vck)[1]
    c = values(vck)[2]
    (nv, nc, nk) = vsize(vck)
    Avck = reshape(A, nv, nc, nk)
    vc_idx = [(iv, ic) for iv in 1:nv, ic in 1:nc]

    norb = numorb(bse)
    AU = Array{ComplexF64}(undef, nk, norb, norb) # ik, ih, ie
    Threads.@threads for I in CartesianIndices(AU)
        (ik, ih, ie) = Tuple(I)
        Uv = conj.(bandk[ik].vectors[ih, v])
        Uc = bandkq[ik].vectors[ie, c]
        AU[I] = sum(vc_idx) do (iv, ic)
            Avck[iv, ic, ik] * Uv[iv] * Uc[ic]
        end
    end

    return RealSpaceBlochExciton_spinblind_wiTDA(q, nk, 1/sqrt(nk), kgrid, AU)
end
