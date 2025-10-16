
struct _cluster_spinless_UJ{
	S <: AbstractReciprocalHoppings,
	L <: AbstractLRCorrection,
	J¹T <: AbstractReciprocalHoppings,
	J²T <: AbstractReciprocalHoppings,
}
	U::UwithLR{S, L}
	J¹::J¹T
	J²::J²T
	UM::Matrix{ComplexF64}
	J¹M::Matrix{ComplexF64}
	J²M::Matrix{ComplexF64}
end

include("./_cluster_spinless_UJ.jl")
include("./_cluster_time_reversal_UJ.jl")

export Kernal_UJ

function Kernal_UJ(
	Kᵈ_U::UwithLR, Kᵈ_J¹::Union{AbstractString, HR}, Kᵈ_J²::Union{AbstractString, HR},
	Kˣ_U::UwithLR, Kˣ_J¹::Union{AbstractString, HR}, Kˣ_J²::Union{AbstractString, HR};
	kwards...)
	if Kᵈ_J¹ isa AbstractString
		Kᵈ_J¹ = ReciprocalHoppings(ReadHR(Kᵈ_J¹))
	elseif J¹ isa HR
		Kᵈ_J¹ = ReciprocalHoppings(Kᵈ_J¹)
	end
	if Kᵈ_J² isa AbstractString
		Kᵈ_J² = ReciprocalHoppings(ReadHR(Kᵈ_J²))
	elseif Kᵈ_J² isa HR
		Kᵈ_J² = ReciprocalHoppings(Kᵈ_J²)
	end
	if Kˣ_J¹ isa AbstractString
		Kˣ_J¹ = ReciprocalHoppings(ReadHR(Kˣ_J¹))
	elseif Kˣ_J¹ isa HR
		Kˣ_J¹ = ReciprocalHoppings(Kˣ_J¹)
	end
	if Kˣ_J² isa AbstractString
		Kˣ_J² = ReciprocalHoppings(ReadHR(Kˣ_J²))
	elseif Kˣ_J² isa HR
		Kˣ_J² = ReciprocalHoppings(Kˣ_J²)
	end
	return Kernal_UJ(Kᵈ_U, Kᵈ_J¹, Kᵈ_J², Kˣ_U, Kˣ_J¹, Kˣ_J²; kwards...)
end
"""

"""
function Kernal_UJ(
	Kᵈ_U::UwithLR{UT_d, CT_d}, Kᵈ_J¹::AbstractReciprocalHoppings{J¹T_d}, Kᵈ_J²::AbstractReciprocalHoppings{J²T_d},
	Kˣ_U::UwithLR{UT_x, CT_x}, Kˣ_J¹::AbstractReciprocalHoppings{J¹T_x}, Kˣ_J²::AbstractReciprocalHoppings{J²T_x};
	upindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
	dnindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
) where {UT_d, CT_d, J¹T_d, J²T_d, UT_x, CT_x, J¹T_x, J²T_x}

	nk = length(kgrid)
	norb = Kᵈ_U.norb

	Kᵈ_UM = Matrix{ComplexF64}(undef, norb, norb)
	Kᵈ_J¹M = Matrix{ComplexF64}(undef, norb, norb)
	Kᵈ_J²M = Matrix{ComplexF64}(undef, norb, norb)
	W = _cluster_spinless_UJ(Kᵈ_U, Kᵈ_J¹, Kᵈ_J², Kᵈ_UM, Kᵈ_J¹M, Kᵈ_J²M)

	Kˣ_UM = Matrix{ComplexF64}(undef, norb, norb)
	Kˣ_J¹M = Matrix{ComplexF64}(undef, norb, norb)
	Kˣ_J²M = Matrix{ComplexF64}(undef, norb, norb)
	V = _cluster_spinless_UJ(Kˣ_U, Kˣ_J¹, Kˣ_J², Kˣ_UM, Kˣ_J¹M, Kˣ_J²M)

	norb_wispin = 2 * norb
	if isnothing(upindex)
		if isnothing(dnindex)
			return _cluster_spinless_UJ_Kernal(W, V)
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
	return _cluster_time_reversal_UJ_Kernal(upindex, dnindex, W, V)
end

function (K::Union{_cluster_spinless_UJ_Kernal, _cluster_time_reversal_UJ_Kernal})(::Val{:initialize})
	K.W.U(K.W.UM, ReducedCoordinates(0,0,0), 1)
	K.W.J¹(K.V.J¹M, ReducedCoordinates(0,0,0))
	K.W.J²(K.V.J²M, ReducedCoordinates(0,0,0))
	K.V.U(K.V.UM, ReducedCoordinates(0,0,0), 1)
	K.V.J¹(K.V.J¹M, ReducedCoordinates(0,0,0))
	K.V.J²(K.V.J²M, ReducedCoordinates(0,0,0))
	return nothing
end
