struct _cluster_time_reversal_UJ_Kernal{
	WT <: _cluster_spinless_UJ,
	VT <: _cluster_spinless_UJ,
} <: AbstractKernalInterAction
	upindex::Vector{Int}
	dnindex::Vector{Int}
	W::WT
	V::VT
end
function (K::_cluster_time_reversal_UJ_Kernal)(ψc′, ψv′, ψv, ψc)::Tuple{ComplexF64, ComplexF64}
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
	U₁ = ψc′_up * transpose(ψc_up) + ψc′_dn * transpose(ψc_dn) #2*n^2
	U₂ = ψv_up * transpose(ψv′_up) + ψv_dn * transpose(ψv′_dn) #2*n^2
	Kᵈ = transpose(diag(U₁)) * K.W.UM * diag(U₂) #n^2+n
	Kᵈ += sum(U₁ .* transpose(U₂) .* K.W.J¹M) #n^2*2
	Kᵈ += sum(U₁ .* U₂ .* K.W.J²M) #n^2*2
	#Kˣ
	U₁ .= ψc′_up * transpose(ψv′_up) + ψc′_dn * transpose(ψv′_dn) #2*n^2
	U₂ .= ψv_up * transpose(ψc_up) + ψv_dn * transpose(ψc_dn) #2*n^2
	Kˣ = transpose(diag(U₁)) * K.V.UM * diag(U₂) #n^2+n
	Kˣ += sum(U₁ .* transpose(U₂) .* K.V.J¹M) #n^2*2
	Kˣ += sum(U₁ .* U₂ .* K.V.J²M) #n^2*2
	#Do not forget this minus.
	return -Kᵈ, Kˣ
end
