struct _spinless_UJA_Kᵈ{
	S <: AbstractReciprocalHoppings,
	L <: AbstractLRCorrection,
	J¹T <: AbstractReciprocalHoppings,
	J²T <: AbstractReciprocalHoppings,
	A¹T <: AbstractReciprocalHoppings,
	A²T <: AbstractReciprocalHoppings,
	A³T <: AbstractReciprocalHoppings,
	A⁴T <: AbstractReciprocalHoppings,
}
	U::UwithLR{S, L}
	J¹::J¹T
	J²::J²T
	A¹::A¹T
	A²::A²T
	A³::A³T
	A⁴::A⁴T
	Uk::Array{ComplexF64, 3}
	J¹q::Matrix{ComplexF64}
	J²k::Array{ComplexF64, 3}
	J²kq::Array{ComplexF64, 3}
	A¹k::Array{ComplexF64, 3}
	A²k::Array{ComplexF64, 3}
	A³k::Array{ComplexF64, 3}
	A⁴k::Array{ComplexF64, 3}
	A¹kq::Array{ComplexF64, 3}
	A⁴kq::Array{ComplexF64, 3}
end
struct _spinless_UJA_Kˣ{
	S <: AbstractReciprocalHoppings,
	L <: AbstractLRCorrection,
	J¹T <: AbstractReciprocalHoppings,
	J²T <: AbstractReciprocalHoppings,
	A¹T <: AbstractReciprocalHoppings,
	A²T <: AbstractReciprocalHoppings,
	A³T <: AbstractReciprocalHoppings,
	A⁴T <: AbstractReciprocalHoppings,
}
	U::UwithLR{S, L}
	J¹::J¹T
	J²::J²T
	A¹::A¹T
	A²::A²T
	A³::A³T
	A⁴::A⁴T
	Uq::Matrix{ComplexF64}
	J¹k::Array{ComplexF64, 3}
	J²k::Array{ComplexF64, 3}
	J²kq::Array{ComplexF64, 3}
	A¹k::Array{ComplexF64, 3}
	A²k::Array{ComplexF64, 3}
	A³k::Array{ComplexF64, 3}
	A⁴k::Array{ComplexF64, 3}
	A¹kq::Array{ComplexF64, 3}
	A³kq::Array{ComplexF64, 3}
end

include("./_spinless_UJA.jl")
include("./_time_reversal_UJA.jl")

export Kernel_UJA

function Kernel_UJA(kgrid::Union{MonkhorstPack, RedKgrid};
	Kᵈ_U::UwithLR,
	Kᵈ_J::Union{AbstractString, HR, AbstractReciprocalHoppings} = "default",
	Kᵈ_J¹::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kᵈ_J,
	Kᵈ_J²::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kᵈ_J,
	Kᵈ_A::Union{AbstractString, HR, AbstractReciprocalHoppings} = "default",
	Kᵈ_A¹::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kᵈ_A,
	Kᵈ_A²::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kᵈ_A,
	Kᵈ_A³::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kᵈ_A,
	Kᵈ_A⁴::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kᵈ_A,
	Kˣ_U::UwithLR,
	Kˣ_J::Union{AbstractString, HR, AbstractReciprocalHoppings} = "default",
	Kˣ_J¹::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kˣ_J,
	Kˣ_J²::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kˣ_J,
	Kˣ_A::Union{AbstractString, HR, AbstractReciprocalHoppings} = "default",
	Kˣ_A¹::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kˣ_A,
	Kˣ_A²::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kˣ_A,
	Kˣ_A³::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kˣ_A,
	Kˣ_A⁴::Union{AbstractString, HR, AbstractReciprocalHoppings} = Kˣ_A,
	upindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
	dnindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
)

	if kgrid isa MonkhorstPack
		kgrid = RedKgrid(kgrid)
	end

	if Kᵈ_J¹ isa AbstractString
		Kᵈ_J¹ = ReciprocalHoppings(ReadHR(Kᵈ_J¹))
	elseif Kᵈ_J¹ isa HR
		Kᵈ_J¹ = ReciprocalHoppings(Kᵈ_J¹)
	end
	if Kᵈ_J² isa AbstractString
		Kᵈ_J² = ReciprocalHoppings(ReadHR(Kᵈ_J²))
	elseif Kᵈ_J² isa HR
		Kᵈ_J² = ReciprocalHoppings(Kᵈ_J²)
	end
	if Kᵈ_A¹ isa AbstractString
		Kᵈ_A¹ = ReciprocalHoppings(ReadHR(Kᵈ_A¹))
	elseif Kᵈ_A¹ isa HR
		Kᵈ_A¹ = ReciprocalHoppings(Kᵈ_A¹)
	end
	if Kᵈ_A² isa AbstractString
		Kᵈ_A² = ReciprocalHoppings(ReadHR(Kᵈ_A²))
	elseif Kᵈ_A² isa HR
		Kᵈ_A² = ReciprocalHoppings(Kᵈ_A²)
	end
	if Kᵈ_A³ isa AbstractString
		Kᵈ_A³ = ReciprocalHoppings(ReadHR(Kᵈ_A³))
	elseif Kᵈ_A³ isa HR
		Kᵈ_A³ = ReciprocalHoppings(Kᵈ_A³)
	end
	if Kᵈ_A⁴ isa AbstractString
		Kᵈ_A⁴ = ReciprocalHoppings(ReadHR(Kᵈ_A⁴))
	elseif Kᵈ_A⁴ isa HR
		Kᵈ_A⁴ = ReciprocalHoppings(Kᵈ_A⁴)
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
	if Kˣ_A¹ isa AbstractString
		Kˣ_A¹ = ReciprocalHoppings(ReadHR(Kˣ_A¹))
	elseif Kˣ_A¹ isa HR
		Kˣ_A¹ = ReciprocalHoppings(Kˣ_A¹)
	end
	if Kˣ_A² isa AbstractString
		Kˣ_A² = ReciprocalHoppings(ReadHR(Kˣ_A²))
	elseif Kˣ_A² isa HR
		Kˣ_A² = ReciprocalHoppings(Kˣ_A²)
	end
	if Kˣ_A³ isa AbstractString
		Kˣ_A³ = ReciprocalHoppings(ReadHR(Kˣ_A³))
	elseif Kˣ_A³ isa HR
		Kˣ_A³ = ReciprocalHoppings(Kˣ_A³)
	end
	if Kˣ_A⁴ isa AbstractString
		Kˣ_A⁴ = ReciprocalHoppings(ReadHR(Kˣ_A⁴))
	elseif Kˣ_A⁴ isa HR
		Kˣ_A⁴ = ReciprocalHoppings(Kˣ_A⁴)
	end

	nk = length(kgrid)
	norb = Kᵈ_U.norb

	Kᵈ_Uk = Array{ComplexF64}(undef, norb, norb, nk)
	Kᵈ_J¹q = Matrix{ComplexF64}(undef, norb, norb)
	Kᵈ_J²k = Array{ComplexF64}(undef, norb, norb, nk)
	Kᵈ_J²kq = Array{ComplexF64}(undef, norb, norb, nk)
	Kᵈ_A¹k = Array{ComplexF64}(undef, norb, norb, nk)
	Kᵈ_A²k = Array{ComplexF64}(undef, norb, norb, nk)
	Kᵈ_A³k = Array{ComplexF64}(undef, norb, norb, nk)
	Kᵈ_A⁴k = Array{ComplexF64}(undef, norb, norb, nk)
	Kᵈ_A¹kq = Array{ComplexF64}(undef, norb, norb, nk)
	Kᵈ_A⁴kq = Array{ComplexF64}(undef, norb, norb, nk)
	W = _spinless_UJA_Kᵈ(Kᵈ_U, Kᵈ_J¹, Kᵈ_J², Kᵈ_A¹, Kᵈ_A², Kᵈ_A³, Kᵈ_A⁴,
		Kᵈ_Uk, Kᵈ_J¹q, Kᵈ_J²k, Kᵈ_J²kq, Kᵈ_A¹k, Kᵈ_A²k, Kᵈ_A³k, Kᵈ_A⁴k, Kᵈ_A¹kq, Kᵈ_A⁴kq)

	Kˣ_Uq = Matrix{ComplexF64}(undef, norb, norb)
	Kˣ_J¹k = Array{ComplexF64}(undef, norb, norb, nk)
	Kˣ_J²k = Array{ComplexF64}(undef, norb, norb, nk)
	Kˣ_J²kq = Array{ComplexF64}(undef, norb, norb, nk)
	Kˣ_A¹k = Array{ComplexF64}(undef, norb, norb, nk)
	Kˣ_A²k = Array{ComplexF64}(undef, norb, norb, nk)
	Kˣ_A³k = Array{ComplexF64}(undef, norb, norb, nk)
	Kˣ_A⁴k = Array{ComplexF64}(undef, norb, norb, nk)
	Kˣ_A¹kq = Array{ComplexF64}(undef, norb, norb, nk)
	Kˣ_A³kq = Array{ComplexF64}(undef, norb, norb, nk)
	V = _spinless_UJA_Kˣ(Kˣ_U, Kˣ_J¹, Kˣ_J², Kˣ_A¹, Kˣ_A², Kˣ_A³, Kˣ_A⁴,
		Kˣ_Uq, Kˣ_J¹k, Kˣ_J²k, Kˣ_J²kq, Kˣ_A¹k, Kˣ_A²k, Kˣ_A³k, Kˣ_A⁴k, Kˣ_A¹kq, Kˣ_A³kq)


	kgrid_Γ = RedKgrid(MonkhorstPack(kgrid.kgrid_size))
	kgrid_addmap = kgridmap(kgrid_Γ.kdirect, kgrid.kdirect, +)
	kgrid_minusmap = kgridmap(kgrid_Γ.kdirect, kgrid.kdirect, -)

	norb_wispin = 2 * norb
	if isnothing(upindex)
		if isnothing(dnindex)
			return _spinless_UJA_Kernel(nk, kgrid, kgrid_Γ, kgrid_addmap, kgrid_minusmap, W, V)
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
	return _time_reversal_UJA_Kernel(nk, kgrid, kgrid_Γ, kgrid_addmap, kgrid_minusmap, upindex, dnindex, W, V)
end


function (K::Union{_spinless_UJA_Kernel, _time_reversal_UJA_Kernel})(::Val{:initialize})
	Threads.@threads for k in Base.OneTo(K.nk)
		Uk = view(K.W.Uk, :, :, k)
		K.W.U(Uk, K.kgrid_Γ[k], K.nk)
		Uk ./= K.nk
		A²k = view(K.W.A²k, :, :, k)
		K.W.A²(A²k, K.kgrid[k], K.nk)
		A²k ./= K.nk
		A³k = view(K.W.A³k, :, :, k)
		K.W.A³(A³k, -K.kgrid[k], K.nk)
		A³k ./= K.nk
		K.V.J¹(K.V.J¹k[k], K.kgrid_Γ[k])
		K.V.J¹k[k] ./= K.nk
		A²k = view(K.V.A²k, :, :, k)
		K.V.A²(A²k, K.kgrid[k], K.nk)
		A²k ./= K.nk
		A⁴k = view(K.V.A⁴k, :, :, k)
		K.V.A⁴(A⁴k, -K.kgrid[k], K.nk)
		A⁴k ./= K.nk
	end
	return nothing
end
function (K::Union{_spinless_UJA_Kernel, _time_reversal_UJA_Kernel})(q::ReducedCoordinates)
	#WJ¹(q)
	K.W.J¹(K.W.J¹q, q)
	K.W.J¹q ./= K.nk
	#VU(q)
	K.V.U(K.V.Uq, q)
	K.V.Uq ./= K.nk
	Threads.@threads for k in Base.OneTo(K.nk)
		#k_Γ+q
		kq = K.kgrid_Γ[k] + q
		#WJ²(k_Γ+q)
		J²kq = view(K.W.J²kq, :, :, k)
		K.W.J²(J²kq, kq)
		J²kq ./= K.nk
		#VJ²(k_Γ+q)
		J²kq = view(K.V.J²kq, :, :, k)
		K.V.J²(J²kq, kq)
		J²kq ./= K.nk
		#k+q
		kq = K.kgrid[k] + q
		A¹kq = view(K.W.A¹kq, :, :, k)
		K.W.A¹(A¹kq, kq)
		A¹kq ./= K.nk
		A⁴kq = view(K.W.A⁴kq, :, :, k)
		K.W.A⁴(A⁴kq, -kq)
		A⁴kq ./= K.nk
		A¹kq = view(K.V.A¹kq, :, :, k)
		K.V.A¹(A¹kq, kq)
		A¹kq ./= K.nk
		A³kq = view(K.V.A³kq, :, :, k)
		K.V.A³(A³kq, -kq)
		A³kq ./= K.nk
	end
	return nothing
end
function (K::Union{_spinless_UJA_Kernel, _time_reversal_UJA_Kernel})(::Val{:initialize_qgrid})
	Threads.@threads for k in Base.OneTo(K.nk)
		k_Γ = K.kgrid_Γ[k]
		J²k = view(K.W.J²k, :, :, k)
		K.W.J²(J²k, k_Γ)
		J²k ./= K.nk
		J²k = view(K.V.J²k, :, :, k)
		K.V.J²(J²k, k_Γ)
		J²k ./= K.nk
		k = K.kgrid[k]
		A¹k = view(K.W.A¹k, :, :, k)
		K.W.A¹(A¹k, k)
		A¹k ./= K.nk
		A⁴k = view(K.W.A⁴k, :, :, k)
		K.W.A⁴(A⁴k, -k)
		A⁴k ./= K.nk
		A¹k = view(K.V.A¹k, :, :, k)
		K.V.A¹(A¹k, k)
		A¹k ./= K.nk
		A³k = view(K.V.A³k, :, :, k)
		K.V.A³(A³k, -k)
		A³k ./= K.nk
	end
	return nothing
end
function (K::Union{_spinless_UJA_Kernel, _time_reversal_UJA_Kernel})(kq_kindex::Vector{Int}, kΓq_kΓindex::Vector{Int}, q::ReducedCoordinates)
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
	#A
	K.W.A¹kq .= K.W.A¹k[:, :, kq_kindex]
	K.W.A⁴kq .= K.W.A⁴k[:, :, kq_kindex]
	K.V.A¹kq .= K.V.A¹k[:, :, kq_kindex]
	K.V.A³kq .= K.V.A³k[:, :, kq_kindex]
	return nothing
end
