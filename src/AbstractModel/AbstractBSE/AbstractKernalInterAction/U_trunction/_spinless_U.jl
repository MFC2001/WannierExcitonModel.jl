
struct _spinless_U_Kernal{
	WT <: _spinless_U_Kᵈ,
	VT <: _spinless_U_Kˣ,
} <: AbstractKernalInterAction
	nk::Int
	kgrid::RedKgrid
	kgrid_Γ::RedKgrid
	kgrid_addmap::Matrix{Int}
	kgrid_minusmap::Matrix{Int}
	W::WT
	V::VT
end
function (K::_spinless_U_Kernal)(k′, k, ψc′, ψv′, ψv, ψc)::Tuple{ComplexF64, ComplexF64}
	conj!(ψc′)
	conj!(ψv)
	U₁ = ψc′ .* ψc #n
	U₂ = ψv .* ψv′ #n
	Kᵈ = transpose(U₁) * K.W.Uk[K.kgrid_minusmap[k′, k]] * U₂ #n^2+n
	U₁ .= ψc′ .* ψv′ #n
	U₂ .= ψv .* ψc #n
	Kˣ = transpose(U₁) * K.V.Uq * U₂ #n^2+n
	#Do not forget this minus.
	return -Kᵈ, Kˣ
end


