struct _time_reversal_UJ_Kernal{
	WT <: _spinless_UJ_Kᵈ,
	VT <: _spinless_UJ_Kˣ,
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
function (K::_time_reversal_UJ_Kernal)(k′, k, ψc′, ψv′, ψv, ψc)::Tuple{ComplexF64, ComplexF64}
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
	Kᵈ = transpose(diag(U₁)) * K.W.Uk[K.kgrid_minusmap[k′, k]] * diag(U₂) #n^2+n
	Kᵈ += sum(U₁ .* transpose(U₂) .* K.W.J¹q) #n^2*2
	Kᵈ += sum(U₁ .* U₂ .* K.W.J²kq[K.kgrid_addmap[k′, k]]) #n^2*2
	#Kˣ
	U₁ .= ψc′_up * transpose(ψv′_up) + ψc′_dn * transpose(ψv′_dn) #2*n^2
	U₂ .= ψv_up * transpose(ψc_up) + ψv_dn * transpose(ψc_dn) #2*n^2
	Kˣ = transpose(diag(U₁)) * K.V.Uq * diag(U₂) #n^2+n
	Kˣ += sum(U₁ .* transpose(U₂) .* K.V.J¹k[K.kgrid_minusmap[k′, k]]) #n^2*2
	Kˣ += sum(U₁ .* U₂ .* K.V.J²kq[K.kgrid_addmap[k′, k]]) #n^2*2
	#Do not forget this minus.
	return -Kᵈ, Kˣ
end
