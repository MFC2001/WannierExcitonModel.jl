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
function realspace(bse::BSEspinblind_wiTDA, q, A::AbstractVector{<:Number})

    kgrid = bse.kgrid
    bandk = bse.bandk
    q = ReducedCoordinates(q)
    bandkq = BAND(map(k -> k + q, kgrid), bse.eleH; vector=true)
    phase_normalize!(:positive_maxabs, bandkq)

    vck = bse.vck
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
