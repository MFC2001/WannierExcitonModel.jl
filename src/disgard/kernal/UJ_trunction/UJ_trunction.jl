struct _spinless_UJ_Kᵈ{
	S <: AbstractReciprocalHoppings,
	L <: AbstractLRCorrection,
	J¹T <: AbstractReciprocalHoppings,
	J²T <: AbstractReciprocalHoppings,
}
	U::UwithLR{S, L}
	J¹::J¹T
	J²::J²T
	Uk::Array{ComplexF64, 3}
	J¹q::Matrix{ComplexF64}
	J²k::Array{ComplexF64, 3}
	J²kq::Array{ComplexF64, 3}
end
struct _spinless_UJ_Kˣ{
	S <: AbstractReciprocalHoppings,
	L <: AbstractLRCorrection,
	J¹T <: AbstractReciprocalHoppings,
	J²T <: AbstractReciprocalHoppings,
}
	U::UwithLR{S, L}
	J¹::J¹T
	J²::J²T
	Uq::Matrix{ComplexF64}
	J¹k::Array{ComplexF64, 3}
	J²k::Array{ComplexF64, 3}
	J²kq::Array{ComplexF64, 3}
end

include("./_spinless_UJ.jl")
include("./_time_reversal_UJ.jl")

export Kernel_UJ

"""
	Kernel_UJ(kgrid::Union{MonkhorstPack, RedKgrid}, 
		Kᵈ_U::UwithLR{UT_d, CT_d}, Kᵈ_J¹::AbstractReciprocalHoppings, Kᵈ_J²::AbstractReciprocalHoppings,
		Kˣ_U::UwithLR{UT_x, CT_x}, Kˣ_J¹::AbstractReciprocalHoppings, Kˣ_J²::AbstractReciprocalHoppings;
		upindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
		dnindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing
	) -> KernelInterAction

Creat a object which can be provided to canstruct a BSE model.
- `kgrid`: should be the same as excitonic basis;
- `K_U`: direct Coulomb interaction with long-range correction between electronic wannier basis;
- `K_J`: exchange Coulomb interaction between electronic wannier basis;
- `up(dn)_index`: the wannier basis index when contains spin.

`Kᵈ` needs screened coulomb interaction between wannier basis, while `Kˣ` needs naked coulomb interaction.
The long-range behavior of direct Coulomb interaction between wannier basis is proportional to \\frac{1}{r}, we can't dismiss it.
Check [`UwithLR`](@ref) to see how to contain the long-range correction.
See more details about `J¹` and `J²` in Theory.

For spinless case, you shouldn't provide `upindex` and `dnindex`.
But for spinful case, you must provide at least one between `upindex` and `dnindex`.
For spinful case, note that the interaction should not contain spin(see more details in Theory.)

	Kernel_UJ(kgrid::Union{MonkhorstPack, RedKgrid}, 
		Kᵈ_U::UwithLR, Kᵈ_J¹::Union{AbstractString, HR}, Kᵈ_J²::Union{AbstractString, HR},
		Kˣ_U::UwithLR, Kˣ_J¹::Union{AbstractString, HR}, Kˣ_J²::Union{AbstractString, HR};
		upindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
		dnindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing
	) -> KernelInterAction
	
This is a shortcut.

!!! note
	For spinful case, we require that the electronic wannier basis satisfy time-reversal symmetry,
	so that we can define which term is `U` and which term is `J`.
"""
function Kernel_UJ(kgrid::Union{MonkhorstPack, RedKgrid};
	Kᵈ_U::UwithLR,
	Kᵈ_J::Union{AbstractString, HR, AbstractReciprocalHoppings} = "default",
	Kᵈ_J¹::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kᵈ_J,
	Kᵈ_J²::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kᵈ_J,
	Kˣ_U::UwithLR,
	Kˣ_J::Union{AbstractString, HR, AbstractReciprocalHoppings} = "default",
	Kˣ_J¹::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kˣ_J,
	Kˣ_J²::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kˣ_J,
	upindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
	dnindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
)

	if kgrid isa MonkhorstPack
		kgrid = RedKgrid(kgrid)
	end

	if Kᵈ_J¹ isa AbstractString
		Kᵈ_J¹ = ReciprocalHoppings(read(Kᵈ_J¹, wannier90_hr))
	elseif Kᵈ_J¹ isa HR
		Kᵈ_J¹ = ReciprocalHoppings(Kᵈ_J¹)
	end
	if Kᵈ_J² isa AbstractString
		Kᵈ_J² = ReciprocalHoppings(read(Kᵈ_J², wannier90_hr))
	elseif Kᵈ_J² isa HR
		Kᵈ_J² = ReciprocalHoppings(Kᵈ_J²)
	end
	if Kˣ_J¹ isa AbstractString
		Kˣ_J¹ = ReciprocalHoppings(read(Kˣ_J¹, wannier90_hr))
	elseif Kˣ_J¹ isa HR
		Kˣ_J¹ = ReciprocalHoppings(Kˣ_J¹)
	end
	if Kˣ_J² isa AbstractString
		Kˣ_J² = ReciprocalHoppings(read(Kˣ_J², wannier90_hr))
	elseif Kˣ_J² isa HR
		Kˣ_J² = ReciprocalHoppings(Kˣ_J²)
	end

	nk = length(kgrid)
	norb = Kᵈ_U.norb

	Kᵈ_Uk = Array{ComplexF64}(undef, norb, norb, nk)
	Kᵈ_J¹q = Matrix{ComplexF64}(undef, norb, norb)
	Kᵈ_J²k = Array{ComplexF64}(undef, norb, norb, nk)
	Kᵈ_J²kq = Array{ComplexF64}(undef, norb, norb, nk)
	W = _spinless_UJ_Kᵈ(Kᵈ_U, Kᵈ_J¹, Kᵈ_J², Kᵈ_Uk, Kᵈ_J¹q, Kᵈ_J²k, Kᵈ_J²kq)

	Kˣ_Uq = Matrix{ComplexF64}(undef, norb, norb)
	Kˣ_J¹k = Array{ComplexF64}(undef, norb, norb, nk)
	Kˣ_J²k = Array{ComplexF64}(undef, norb, norb, nk)
	Kˣ_J²kq = Array{ComplexF64}(undef, norb, norb, nk)
	V = _spinless_UJ_Kˣ(Kˣ_U, Kˣ_J¹, Kˣ_J², Kˣ_Uq, Kˣ_J¹k, Kˣ_J²k, Kˣ_J²kq)


	kgrid_Γ = RedKgrid(MonkhorstPack(kgrid.kgrid_size))
	kgrid_addmap = kgridmap(kgrid_Γ.kdirect, kgrid.kdirect, +)
	kgrid_minusmap = kgridmap(kgrid_Γ.kdirect, kgrid.kdirect, -)

	norb_wispin = 2 * norb
	if isnothing(upindex)
		if isnothing(dnindex)
			return _spinless_UJ_Kernel(nk, kgrid, kgrid_Γ, kgrid_addmap, kgrid_minusmap, W, V)
		else
			dnindex = collect(Int.(dnindex))
			upindex = setdiff(1:norb_wispin, dnindex)
		end
	else
		upindex = collect(Int.(upindex))
		if isnothing(dnindex)
			dnindex = setdiff(1:norb_wispin, upindex)
		else
			dnindex = collect(Int.(dnindex))
		end
	end
	return _time_reversal_UJ_Kernel(nk, kgrid, kgrid_Γ, kgrid_addmap, kgrid_minusmap, upindex, dnindex, W, V)
end

function (K::Union{_spinless_UJ_Kernel, _time_reversal_UJ_Kernel})(::Val{:initialize})
	Threads.@threads for k in Base.OneTo(K.nk)
		Uk = view(K.W.Uk, :, :, k)
		K.W.U(Uk, K.kgrid_Γ[k], K.nk)
		Uk ./= K.nk
		J¹k = view(K.V.J¹k, :, :, k)
		K.V.J¹(J¹k, K.kgrid_Γ[k])
		J¹k ./= K.nk
	end
	return nothing
end
function (K::Union{_spinless_UJ_Kernel, _time_reversal_UJ_Kernel})(q::ReducedCoordinates)
	#WJ¹(q)
	K.W.J¹(K.W.J¹q, q)
	K.W.J¹q ./= K.nk
	#VU(q)
	K.V.U(K.V.Uq, q)
	K.V.Uq ./= K.nk
	Threads.@threads for k in Base.OneTo(K.nk)
		kq = K.kgrid_Γ[k] + q
		#WJ²(k_Γ+q)
		J²kq = view(K.W.J²kq, :, :, k)
		K.W.J²(J²kq, kq)
		J²kq ./= K.nk
		#VJ²(k_Γ+q)
		J²kq = view(K.V.J²kq, :, :, k)
		K.V.J²(J²kq, kq)
		J²kq ./= K.nk
	end
	return nothing
end
function (K::Union{_spinless_UJ_Kernel, _time_reversal_UJ_Kernel})(::Val{:initialize_qgrid})
	Threads.@threads for k in Base.OneTo(K.nk)
		J²k = view(K.W.J²k, :, :, k)
		K.W.J²(J²k, K.kgrid_Γ[k])
		J²k ./= K.nk
		J²k = view(K.V.J²k, :, :, k)
		K.V.J²(J²k, K.kgrid_Γ[k])
		J²k ./= K.nk
	end
	return nothing
end
function (K::Union{_spinless_UJ_Kernel, _time_reversal_UJ_Kernel})(kq_kindex::Vector{Int}, kΓq_kΓindex::Vector{Int}, q::ReducedCoordinates)
	#WJ¹(q)
	K.W.J¹(K.J¹q, q)
	K.J¹q ./= K.nk
	#WJ²(k_Γ+q)
	K.W.J²kq .= K.W.J²k[:, :, kΓq_kΓindex]
	#VU(q)
	K.V.U(K.V.Uq, q)
	K.V.Uq ./= K.nk
	#VJ²(k_Γ+q)
	K.V.J²kq .= K.V.J²k[:, :, kΓq_kΓindex]
	return nothing
end
