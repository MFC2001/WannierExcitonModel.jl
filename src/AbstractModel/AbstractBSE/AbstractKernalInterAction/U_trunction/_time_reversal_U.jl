struct _time_reversal_U_Kernal{
	WT <: _spinless_U_Kᵈ,
	VT <: _spinless_U_Kˣ,
} <: AbstractKernalInterAction
	nk::Int
	kgrid::RedKgrid
	kgrid_Γ::RedKgrid
	kgrid_addmap::Matrix{Int}
	kgrid_minusmap::Matrix{Int}
	upindex::Vector{Int}
	dnindex::Vector{Int}
	W::WT
	V::VT
end
function (K::_time_reversal_U_Kernal)(k′, k, ψc′, ψv′, ψv, ψc)::Tuple{ComplexF64, ComplexF64}
	conj!(ψc′)
	conj!(ψv)
	@views begin
		ψc′_up = ψc′[K.upindex]
		ψc′_dn = ψc′[K.dnindex]
		ψv′_up = ψv′[K.upindex]
		ψv′_dn = ψv′[K.dnindex]
		ψv_up = ψv[K.upindex]
		ψv_dn = ψv[K.dnindex]
		ψc_up = ψc[K.upindex]
		ψc_dn = ψc[K.dnindex]
	end
	#Kᵈ
	U₁ = ψc′_up .* ψc_up + ψc′_dn .* ψc_dn
	U₂ = ψv_up .* ψv′_up + ψv_dn .* ψv′_dn
	Kᵈ = transpose(U₁) * K.W.Uk[K.kgrid_minusmap[k′, k]] * U₂
	#Kˣ
	U₁ .= ψc′_up .* ψv′_up + ψc′_dn .* ψv′_dn
	U₂ .= ψv_up .* ψc_up + ψv_dn .* ψc_dn
	Kˣ = transpose(U₁) * K.V.Uq * U₂
	#Do not forget this minus.
	return -Kᵈ, Kˣ
end
