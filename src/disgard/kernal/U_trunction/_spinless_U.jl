
struct _spinless_U_Kernel{
	WT <: _spinless_U_Kᵈ,
	VT <: _spinless_U_Kˣ,
} <: KernelInterAction
	nk::Int
	kgrid::RedKgrid
	kgrid_Γ::RedKgrid
	kgrid_addmap::Matrix{Int}
	kgrid_minusmap::Matrix{Int}
	W::WT
	V::VT
end
function (K::_spinless_U_Kernel)(k′, k, ψc′, ψv′, ψv, ψc)::Tuple{ComplexF64, ComplexF64}
	conj!(ψc′)
	conj!(ψv)
	U₁ = ψc′ .* ψc #n
	U₂ = ψv .* ψv′ #n
	Ukk = view(K.W.Uk, :, :, K.kgrid_minusmap[k′, k])
	Kᵈ = transpose(U₁) * Ukk * U₂ #n^2+n
	U₁ .= ψc′ .* ψv′ #n
	U₂ .= ψv .* ψc #n
	Kˣ = transpose(U₁) * K.V.Uq * U₂ #n^2+n
	#Do not forget this minus.
	return -Kᵈ, Kˣ
end


