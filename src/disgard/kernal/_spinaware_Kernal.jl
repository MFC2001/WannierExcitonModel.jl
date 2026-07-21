struct _spinaware_Kernel{
	WT <: Union{_U_Kᵈ, _U_pair_Kᵈ},
	VT <: Union{_U_Kˣ, _U_pair_Kˣ},
} <: KernelInterAction
	nk::Int
	kgrid::RedKgrid
	upindex::Vector{Int}
	dnindex::Vector{Int}
	W_pairupindex::Vector{Int}
	W_pairdnindex::Vector{Int}
	V_pairupindex::Vector{Int}
	V_pairdnindex::Vector{Int}
	W::WT
	V::VT
end
function (K::_spinaware_Kernel)(ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64}
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
	Kᵈ = K.W(ψv′_up, ψc′_up, k′, ψv_up, ψc_up, k)
	Kᵈ += K.W(ψv′_up, ψc′_dn, k′, ψv_up, ψc_dn, k)
	Kᵈ += K.W(ψv′_dn, ψc′_up, k′, ψv_dn, ψc_up, k)
	Kᵈ += K.W(ψv′_dn, ψc′_dn, k′, ψv_dn, ψc_dn, k)
	Kˣ = K.V(ψv′_up, ψc′_up, k′, ψv_up, ψc_up, k)
	Kˣ += K.V(ψv′_up, ψc′_up, k′, ψv_dn, ψc_dn, k)
	Kˣ += K.V(ψv′_dn, ψc′_dn, k′, ψv_up, ψc_up, k)
	Kˣ += K.V(ψv′_dn, ψc′_dn, k′, ψv_dn, ψc_dn, k)
	return Kᵈ, Kˣ
end
function (K::_spinaware_Kernel)(buffer_Kᵈ, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64}
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
	Kᵈ = K.W(buffer_Kᵈ..., ψv′_up, ψc′_up, k′, ψv_up, ψc_up, k)
	Kᵈ += K.W(buffer_Kᵈ..., ψv′_up, ψc′_dn, k′, ψv_up, ψc_dn, k)
	Kᵈ += K.W(buffer_Kᵈ..., ψv′_dn, ψc′_up, k′, ψv_dn, ψc_up, k)
	Kᵈ += K.W(buffer_Kᵈ..., ψv′_dn, ψc′_dn, k′, ψv_dn, ψc_dn, k)
	Kˣ = K.V(buffer_Kˣ..., ψv′_up, ψc′_up, k′, ψv_up, ψc_up, k)
	Kˣ += K.V(buffer_Kˣ..., ψv′_up, ψc′_up, k′, ψv_dn, ψc_dn, k)
	Kˣ += K.V(buffer_Kˣ..., ψv′_dn, ψc′_dn, k′, ψv_up, ψc_up, k)
	Kˣ += K.V(buffer_Kˣ..., ψv′_dn, ψc′_dn, k′, ψv_dn, ψc_dn, k)
	return Kᵈ, Kˣ
end
# work for Sz-conserved BSE, when interaction is spinful.
function (K::_spinaware_Kernel)(::Val{:uu}, ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64}
	Kᵈ = K.W(ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.upindex, iw2 = K.upindex, ip1 = K.W_pairupindex, ip2 = K.W_pairupindex)
	Kˣ = K.V(ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.upindex, iw2 = K.upindex, ip1 = K.V_pairupindex, ip2 = K.V_pairupindex)
	return Kᵈ, Kˣ
end
function (K::_spinaware_Kernel)(::Val{:uu}, buffer_Kᵈ, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64}
	Kᵈ = K.W(buffer_Kᵈ..., ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.upindex, iw2 = K.upindex, ip1 = K.W_pairupindex, ip2 = K.W_pairupindex)
	Kˣ = K.V(buffer_Kˣ..., ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.upindex, iw2 = K.upindex, ip1 = K.V_pairupindex, ip2 = K.V_pairupindex)
	return Kᵈ, Kˣ
end
function (K::_spinaware_Kernel)(::Val{:dd}, ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64}
	Kᵈ = K.W(ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.dnindex, iw2 = K.dnindex, ip1 = K.W_pairdnindex, ip2 = K.W_pairdnindex)
	Kˣ = K.V(ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.dnindex, iw2 = K.dnindex, ip1 = K.V_pairdnindex, ip2 = K.V_pairdnindex)
	return Kᵈ, Kˣ
end
function (K::_spinaware_Kernel)(::Val{:dd}, buffer_Kᵈ, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64}
	Kᵈ = K.W(buffer_Kᵈ..., ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.dnindex, iw2 = K.dnindex, ip1 = K.W_pairdnindex, ip2 = K.W_pairdnindex)
	Kˣ = K.V(buffer_Kˣ..., ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.dnindex, iw2 = K.dnindex, ip1 = K.V_pairdnindex, ip2 = K.V_pairdnindex)
	return Kᵈ, Kˣ
end
function (K::_spinaware_Kernel)(::Val{:uudd}, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64
	Kˣ = K.V(ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.upindex, iw2 = K.dnindex, ip1 = K.V_pairupindex, ip2 = K.V_pairdnindex)
	return Kˣ
end
function (K::_spinaware_Kernel)(::Val{:uudd}, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64
	Kˣ = K.V(buffer_Kˣ..., ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.upindex, iw2 = K.dnindex, ip1 = K.V_pairupindex, ip2 = K.V_pairdnindex)
	return Kˣ
end
function (K::_spinaware_Kernel)(::Val{:dduu}, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64
	Kˣ = K.V(ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.dnindex, iw2 = K.upindex, ip1 = K.V_pairdnindex, ip2 = K.V_pairupindex)
	return Kˣ
end
function (K::_spinaware_Kernel)(::Val{:dduu}, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64
	Kˣ = K.V(buffer_Kˣ..., ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.dnindex, iw2 = K.upindex, ip1 = K.V_pairdnindex, ip2 = K.V_pairupindex)
	return Kˣ
end
function (K::_spinaware_Kernel)(::Val{:ud}, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64
	Kᵈ = K.W(ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.dnindex, iw2 = K.upindex, ip1 = K.W_pairdnindex, ip2 = K.W_pairupindex)
	return Kᵈ
end
function (K::_spinaware_Kernel)(::Val{:ud}, buffer_Kᵈ, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64
	Kᵈ = K.W(buffer_Kᵈ..., ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.dnindex, iw2 = K.upindex, ip1 = K.W_pairdnindex, ip2 = K.W_pairupindex)
	return Kᵈ
end
function (K::_spinaware_Kernel)(::Val{:du}, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64
	Kᵈ = K.W(ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.upindex, iw2 = K.dnindex, ip1 = K.W_pairupindex, ip2 = K.W_pairdnindex)
	return Kᵈ
end
function (K::_spinaware_Kernel)(::Val{:du}, buffer_Kᵈ, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64
	Kᵈ = K.W(buffer_Kᵈ..., ψv′, ψc′, k′, ψv, ψc, k; iw1 = K.upindex, iw2 = K.dnindex, ip1 = K.W_pairupindex, ip2 = K.W_pairdnindex)
	return Kᵈ
end
