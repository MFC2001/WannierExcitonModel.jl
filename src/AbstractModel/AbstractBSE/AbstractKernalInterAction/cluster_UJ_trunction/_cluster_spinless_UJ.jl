
struct _cluster_spinless_UJ_Kernal{
	WT <: _cluster_spinless_UJ,
	VT <: _cluster_spinless_UJ,
} <: AbstractKernalInterAction
	W::WT
	V::VT
end
function (K::_cluster_spinless_UJ_Kernal)(ψc′, ψv′, ψv, ψc)::Tuple{ComplexF64, ComplexF64}
	conj!(ψc′)
	conj!(ψv)
	U₁ = ψc′ .* ψc #n
	U₂ = ψv .* ψv′ #n
	Kᵈ = transpose(U₁) * K.W.UM[K.kgrid_minusmap[k′, k]] * U₂ #n^2+n
	Kˣ = transpose(U₁) * K.V.J¹M[K.kgrid_minusmap[k′, k]] * U₂ #n^2+n
	U₁ .= ψc′ .* ψv′ #n
	U₂ .= ψv .* ψc #n
	Kᵈ += transpose(U₁) * K.W.J¹M * U₂ #n^2+n
	Kˣ += transpose(U₁) * K.V.UM * U₂ #n^2+n
	U₁ .= ψc′ .* ψv #n
	U₂ .= ψv′ .* ψc #n
	Kᵈ += transpose(U₁) * K.W.J²M * U₂ #n^2+n
	Kˣ += transpose(U₁) * K.V.J²M * U₂ #n^2+n
	#Do not forget this minus.
	return -Kᵈ, Kˣ
end


