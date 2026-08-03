struct _spinaware_Kernel{
	WT <: Union{_spinblind_Kᵈ, _spinaware_Kᵈ},
	VT <: Union{_spinblind_Kˣ, _spinaware_Kˣ},
} <: KernelInterAction
	nk::Int
	kgrid::RedKgrid
	upindex::Vector{Int}
	dnindex::Vector{Int}
	W::WT
	V::VT
end
function (K::_spinaware_Kernel{WT, VT})(ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64} where {WT <: _spinblind_Kᵈ, VT <: _spinblind_Kˣ}
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
function (K::_spinaware_Kernel{WT, VT})(buffer_Kᵈ, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64} where {WT <: _spinblind_Kᵈ, VT <: _spinblind_Kˣ}
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
function (K::_spinaware_Kernel{WT, VT})(::Val{:uu}, ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64} where {WT <: _spinaware_Kᵈ, VT <: _spinaware_Kˣ}
	Kᵈ = K.W(Val(:uu), ψv′, ψc′, k′, ψv, ψc, k)
	Kˣ = K.V(Val(:uu), ψv′, ψc′, k′, ψv, ψc, k)
	return Kᵈ, Kˣ
end
function (K::_spinaware_Kernel{WT, VT})(::Val{:uu}, buffer_Kᵈ, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64} where {WT <: _spinaware_Kᵈ, VT <: _spinaware_Kˣ}
	Kᵈ = K.W(Val(:uu), buffer_Kᵈ..., ψv′, ψc′, k′, ψv, ψc, k)
	Kˣ = K.V(Val(:uu), buffer_Kˣ..., ψv′, ψc′, k′, ψv, ψc, k)
	return Kᵈ, Kˣ
end
function (K::_spinaware_Kernel{WT, VT})(::Val{:dd}, ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64} where {WT <: _spinaware_Kᵈ, VT <: _spinaware_Kˣ}
	Kᵈ = K.W(Val(:dd), ψv′, ψc′, k′, ψv, ψc, k)
	Kˣ = K.V(Val(:dd), ψv′, ψc′, k′, ψv, ψc, k)
	return Kᵈ, Kˣ
end
function (K::_spinaware_Kernel{WT, VT})(::Val{:dd}, buffer_Kᵈ, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)::Tuple{ComplexF64, ComplexF64} where {WT <: _spinaware_Kᵈ, VT <: _spinaware_Kˣ}
	Kᵈ = K.W(Val(:dd), buffer_Kᵈ..., ψv′, ψc′, k′, ψv, ψc, k)
	Kˣ = K.V(Val(:dd), buffer_Kˣ..., ψv′, ψc′, k′, ψv, ψc, k)
	return Kᵈ, Kˣ
end
function (K::_spinaware_Kernel{WT, VT})(::Val{:uudd}, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64 where {WT <: _spinaware_Kᵈ, VT <: _spinaware_Kˣ}
	Kˣ = K.V(Val(:uudd), ψv′, ψc′, k′, ψv, ψc, k)
	return Kˣ
end
function (K::_spinaware_Kernel{WT, VT})(::Val{:uudd}, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64 where {WT <: _spinaware_Kᵈ, VT <: _spinaware_Kˣ}
	Kˣ = K.V(Val(:uudd), buffer_Kˣ..., ψv′, ψc′, k′, ψv, ψc, k)
	return Kˣ
end
function (K::_spinaware_Kernel{WT, VT})(::Val{:dduu}, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64 where {WT <: _spinaware_Kᵈ, VT <: _spinaware_Kˣ}
	Kˣ = K.V(Val(:dduu), ψv′, ψc′, k′, ψv, ψc, k)
	return Kˣ
end
function (K::_spinaware_Kernel{WT, VT})(::Val{:dduu}, buffer_Kˣ, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64 where {WT <: _spinaware_Kᵈ, VT <: _spinaware_Kˣ}
	Kˣ = K.V(Val(:dduu), buffer_Kˣ..., ψv′, ψc′, k′, ψv, ψc, k)
	return Kˣ
end
function (K::_spinaware_Kernel{WT, VT})(::Val{:ud}, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64 where {WT <: _spinaware_Kᵈ, VT <: _spinaware_Kˣ}
	Kᵈ = K.W(Val(:ud), ψv′, ψc′, k′, ψv, ψc, k)
	return Kᵈ
end
function (K::_spinaware_Kernel{WT, VT})(::Val{:ud}, buffer_Kᵈ, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64 where {WT <: _spinaware_Kᵈ, VT <: _spinaware_Kˣ}
	Kᵈ = K.W(Val(:ud), buffer_Kᵈ..., ψv′, ψc′, k′, ψv, ψc, k)
	return Kᵈ
end
function (K::_spinaware_Kernel{WT, VT})(::Val{:du}, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64 where {WT <: _spinaware_Kᵈ, VT <: _spinaware_Kˣ}
	Kᵈ = K.W(Val(:du), ψv′, ψc′, k′, ψv, ψc, k)
	return Kᵈ
end
function (K::_spinaware_Kernel{WT, VT})(::Val{:du}, buffer_Kᵈ, ψv′, ψc′, k′, ψv, ψc, k)::ComplexF64 where {WT <: _spinaware_Kᵈ, VT <: _spinaware_Kˣ}
	Kᵈ = K.W(Val(:du), buffer_Kᵈ..., ψv′, ψc′, k′, ψv, ψc, k)
	return Kᵈ
end
