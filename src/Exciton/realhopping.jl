export Kᵈhopping, Kˣhopping
function Kᵈhopping(Kᵈ_U, Kᵈ_pair, p₁, p₂, R₀)
    Kᵈ_p₁ = (p₁[2], p₂[2], R₀ + p₂[3] - p₁[3])
    Kᵈ_p₂ = (p₁[1], p₂[1], R₀)
    if WannierInterActionBase._isonsite_pair(Kᵈ_p₁) && WannierInterActionBase._isonsite_pair(Kᵈ_p₂)
        return - Kᵈ_U[p₁[2], p₁[1], -p₁[3]]
    else
        return - Kᵈ_pair[Kᵈ_p₁, Kᵈ_p₂, -p₁[3]]
    end
end
function Kˣhopping(Kˣ_U, Kˣ_pair, p₁, p₂, R₀)
    if WannierInterActionBase._isonsite_pair(p₁) && WannierInterActionBase._isonsite_pair(p₂)
        return conj(Kˣ_U[p₁[1], p₂[1], R₀])
    else
        return conj(Kˣ_pair[p₁, p₂, R₀])
    end
end
