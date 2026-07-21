struct _time_reversal_UJA_Kernel{
	WT <: _spinless_UJ_Kᵈ,
	VT <: _spinless_UJ_Kˣ,
} <: KernelInterAction
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
function (K::_time_reversal_UJA_Kernel)(k′, k, ψc′, ψv′, ψv, ψc)::Tuple{ComplexF64, ComplexF64}
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
	U₁_diag = diag(U₁)
	U₂_diag = diag(U₂)
	mat_buffer = similar(U₁, ComplexF64)
	vec_buffer = Vector{ComplexF64}(undef, size(mat_buffer, 1))
	Ukk = view(K.W.Uk, :, :, K.kgrid_minusmap[k′, k])
	mul!(vec_buffer, Ukk, U₂_diag)
	vec_buffer .*= U₁_diag
	Kᵈ = sum(vec_buffer) #n^2+n
	mat_buffer .= U₁ .* transpose(U₂) .* K.W.J¹q
	Kᵈ += sum(mat_buffer) #n^2*2
	J²kkq = view(K.W.J²kq, :, :, K.kgrid_addmap[k′, k])
	mat_buffer .= U₁ .* U₂ .* J²kkq
	Kᵈ += sum(mat_buffer) #n^2*2
	#A
	A¹kq = view(K.W.A¹kq, :, :, k′)
	mat_buffer .= U₁ .* A¹kq
	mul!(vec_buffer, mat_buffer, U₂_diag)
	Kᵈ += sum(vec_buffer) #n^2*2
	A²k = view(K.W.A²k, :, :, k)
	mat_buffer .= U₂ .* A²k
	mul!(vec_buffer, mat_buffer, U₁_diag)
	Kᵈ += sum(vec_buffer) #n^2*2
	A³k = view(K.W.A³k, :, :, k′)
	mat_buffer .= transpose(U₂) .* A³k
	mul!(vec_buffer, mat_buffer, U₁_diag)
	Kᵈ += sum(vec_buffer) #n^2*2
	A⁴kq = view(K.W.A⁴kq, :, :, k)
	mat_buffer .= transpose(U₁) .* A⁴kq
	mul!(vec_buffer, mat_buffer, U₂_diag)
	Kᵈ += sum(vec_buffer) #n^2*2
	#Kˣ
	U₁ .= ψc′_up * transpose(ψv′_up) + ψc′_dn * transpose(ψv′_dn) #2*n^2
	U₂ .= ψv_up * transpose(ψc_up) + ψv_dn * transpose(ψc_dn) #2*n^2
	U₁_diag = diag(U₁)
	U₂_diag = diag(U₂)
	mul!(vec_buffer, K.V.Uq, U₂_diag)
	vec_buffer .*= U₁_diag
	Kˣ = sum(vec_buffer) #n^2+n
	J¹kk = view(K.V.J¹k, :, :, K.kgrid_minusmap[k′, k])
	mat_buffer .= U₁ .* transpose(U₂) .* J¹kk
	Kˣ += sum(mat_buffer) #n^2*2
	J²kkq = view(K.V.J²kq, :, :, K.kgrid_addmap[k′, k])
	mat_buffer .= U₁ .* U₂ .* J²kkq
	Kˣ += sum(mat_buffer) #n^2*2
	#A
	A¹kq = view(K.V.A¹kq, :, :, k′)
	mat_buffer .= U₁ .* A¹kq
	mul!(vec_buffer, mat_buffer, U₂_diag)
	Kˣ += sum(vec_buffer) #n^2*2
	A²k = view(K.V.A²k, :, :, k)
	mat_buffer .= U₂ .* A²k
	mul!(vec_buffer, mat_buffer, U₁_diag)
	Kˣ += sum(vec_buffer) #n^2*2
	A³kq = view(K.V.A³kq, :, :, k)
	mat_buffer .= transpose(U₂) .* A³kq
	mul!(vec_buffer, mat_buffer, U₁_diag)
	Kˣ += sum(vec_buffer) #n^2*2
	A⁴k = view(K.V.A⁴k, :, :, k′)
	mat_buffer .= transpose(U₁) .* A⁴k
	mul!(vec_buffer, mat_buffer, U₂_diag)
	Kˣ += sum(vec_buffer) #n^2*2
	#Do not forget this minus.
	return -Kᵈ, Kˣ
end
