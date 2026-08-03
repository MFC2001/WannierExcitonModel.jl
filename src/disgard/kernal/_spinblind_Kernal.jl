struct _spinblind_Kernel{
	WT <: Union{_U_Kᵈ, _U_pair_Kᵈ},
	VT <: Union{_U_Kˣ, _U_pair_Kˣ},
} <: KernelInterAction
	nk::Int
	kgrid::RedKgrid
	W::WT
	V::VT
end
function (K::_spinblind_Kernel)(ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64}
	Kᵈ = K.W(ψv′, ψc′, k′, ψv, ψc, k)
	Kˣ = K.V(ψv′, ψc′, k′, ψv, ψc, k)
	return Kᵈ, Kˣ
end
function (K::_spinblind_Kernel)(buffer_Kᵈ, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64}
	Kᵈ = K.W(buffer_Kᵈ..., ψv′, ψc′, k′, ψv, ψc, k)
	Kˣ = K.V(buffer_Kˣ..., ψv′, ψc′, k′, ψv, ψc, k)
	return Kᵈ, Kˣ
end
# work for Sz-conserved BSE, when interaction is spinless.
function (K::_spinblind_Kernel)(::Union{Val{:uu}, Val{:dd}}, ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64}
	Kᵈ = K.W(ψv′, ψc′, k′, ψv, ψc, k)
	Kˣ = K.V(ψv′, ψc′, k′, ψv, ψc, k)
	return Kᵈ, Kˣ
end
function (K::_spinblind_Kernel)(::Union{Val{:uu}, Val{:dd}}, buffer_Kᵈ, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64}
	Kᵈ = K.W(buffer_Kᵈ..., ψv′, ψc′, k′, ψv, ψc, k)
	Kˣ = K.V(buffer_Kˣ..., ψv′, ψc′, k′, ψv, ψc, k)
	return Kᵈ, Kˣ
end
function (K::_spinblind_Kernel)(::Union{Val{:uudd}, Val{:dduu}}, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64
	Kˣ = K.V(ψv′, ψc′, k′, ψv, ψc, k)
	return Kˣ
end
function (K::_spinblind_Kernel)(::Union{Val{:uudd}, Val{:dduu}}, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64
	Kˣ = K.V(buffer_Kˣ..., ψv′, ψc′, k′, ψv, ψc, k)
	return Kˣ
end
function (K::_spinblind_Kernel)(::Union{Val{:ud}, Val{:du}}, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64
	Kᵈ = K.W(ψv′, ψc′, k′, ψv, ψc, k)
	return Kᵈ
end
function (K::_spinblind_Kernel)(::Union{Val{:ud}, Val{:du}}, buffer_Kᵈ, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64
	Kᵈ = K.W(buffer_Kᵈ..., ψv′, ψc′, k′, ψv, ψc, k)
	return Kᵈ
end
