export Kernel
"""
	abstract type KernelInterAction <: AbstractInterAction end

This abstract type is used to calculate kernel term in excitonic BSE model.

```julia
>julia 

```
"""
abstract type KernelInterAction <: AbstractInterAction end

function (K::KernelInterAction)(sym::Symbol, args...; kwargs...)
	return K(Val(sym), args...; kwargs...)
end
(K::KernelInterAction)(q::AbstractVector{<:Real}; isΓ::Bool = false) =
	K(ReducedCoordinates(q); isΓ)
function (K::KernelInterAction)(q::ReducedCoordinates; isΓ::Bool = false)
	K.W(q)
	K.V(q; isΓ)
	return nothing
end
function (K::KernelInterAction)(::Val{:buffer})
	buffer_Kᵈ = K.W(Val(:buffer))
	buffer_Kˣ = K.V(Val(:buffer))
	return buffer_Kᵈ, buffer_Kˣ
end
(K::KernelInterAction)(::Val{:buffer_Kᵈ}) = K.W(Val(:buffer))
(K::KernelInterAction)(::Val{:buffer_Kˣ}) = K.V(Val(:buffer))
function (K::KernelInterAction)(::Val{:buffer_uu})
	buffer_Kᵈ = K.W(Val(:buffer_uu))
	buffer_Kˣ = K.V(Val(:buffer_uu))
	return buffer_Kᵈ, buffer_Kˣ
end
function (K::KernelInterAction)(::Val{:buffer_dd})
	buffer_Kᵈ = K.W(Val(:buffer_dd))
	buffer_Kˣ = K.V(Val(:buffer_dd))
	return buffer_Kᵈ, buffer_Kˣ
end
function (K::KernelInterAction)(::Val{:buffer_uudd})
	buffer_Kˣ = K.V(Val(:buffer_uudd))
	return buffer_Kˣ
end
function (K::KernelInterAction)(::Val{:buffer_dduu})
	buffer_Kˣ = K.V(Val(:buffer_dduu))
	return buffer_Kˣ
end
function (K::KernelInterAction)(::Val{:buffer_ud})
	buffer_Kᵈ = K.W(Val(:buffer_ud))
	return buffer_Kᵈ
end
function (K::KernelInterAction)(::Val{:buffer_du})
	buffer_Kᵈ = K.W(Val(:buffer_du))
	return buffer_Kᵈ
end

include("./spinblind/spinblind.jl")
include("./spinaware/spinaware.jl")
include("./Kernel/Kernel.jl")

function Kernel(kgrid::Union{MonkhorstPack, RedKgrid},
	Kᵈ_U::DirectInterAction,
	Kˣ_U::DirectInterAction;
	interaction::Symbol = :spinblind, atol::Union{Real, Nothing} = nothing,
	Kᵈ_J::Union{HR, AbstractReciprocalHoppings, Nothing} = nothing,
	Kᵈ_J¹::Union{HR, AbstractReciprocalHoppings, Nothing} = Kᵈ_J,
	Kᵈ_J²::Union{HR, AbstractReciprocalHoppings, Nothing} = Kᵈ_J,
	Kᵈ_pair::Union{WanIntijklR, Nothing} = nothing,
	Kˣ_J::Union{HR, AbstractReciprocalHoppings, Nothing} = nothing,
	Kˣ_J¹::Union{HR, AbstractReciprocalHoppings, Nothing} = Kˣ_J,
	Kˣ_J²::Union{HR, AbstractReciprocalHoppings, Nothing} = Kˣ_J,
	Kˣ_pair::Union{WanIntijklR, Nothing} = nothing,
	upindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
	dnindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
)

	if kgrid isa MonkhorstPack
		kgrid = RedKgrid(kgrid)
	end
	if !isnothing(atol)
		if iszero(atol)
			atol = nothing
		else
			@assert atol > 0 "atol should be positive!"
		end
	end
	if Kᵈ_J¹ isa HR
		Kᵈ_J¹ = ReciprocalHoppings(Kᵈ_J¹)
	end
	if Kᵈ_J¹ isa AbstractReciprocalHoppings && atol isa Number
		Kᵈ_J¹ = prune(Kᵈ_J¹, atol)
	end
	if Kᵈ_J² isa HR
		Kᵈ_J² = ReciprocalHoppings(Kᵈ_J²)
	end
	if Kᵈ_J² isa AbstractReciprocalHoppings && atol isa Number
		Kᵈ_J² = prune(Kᵈ_J², atol)
	end
	if Kᵈ_pair isa WanIntijklR && atol isa Number
		Kᵈ_pair = prune(Kᵈ_pair, atol)
	end
	if Kˣ_J¹ isa HR
		Kˣ_J¹ = ReciprocalHoppings(Kˣ_J¹)
	end
	if Kˣ_J¹ isa AbstractReciprocalHoppings && atol isa Number
		Kˣ_J¹ = prune(Kˣ_J¹, atol)
	end
	if Kˣ_J² isa HR
		Kˣ_J² = ReciprocalHoppings(Kˣ_J²)
	end
	if Kˣ_J² isa AbstractReciprocalHoppings && atol isa Number
		Kˣ_J² = prune(Kˣ_J², atol)
	end
	if Kˣ_pair isa WanIntijklR && atol isa Number
		Kˣ_pair = prune(Kˣ_pair, atol)
	end

	return Kernel(Val(interaction), kgrid, Kᵈ_U, Kˣ_U;
		Kᵈ_J¹, Kᵈ_J², Kᵈ_pair, Kˣ_J¹, Kˣ_J², Kˣ_pair,
		upindex, dnindex)
end

function Kernel(::Val{:spinblind}, kgrid::RedKgrid,
	Kᵈ_U::DirectInterAction,
	Kˣ_U::DirectInterAction;
	Kᵈ_J¹::Union{AbstractReciprocalHoppings, Nothing} = nothing,
	Kᵈ_J²::Union{AbstractReciprocalHoppings, Nothing} = nothing,
	Kᵈ_pair::Union{WanIntijklR, Nothing} = nothing,
	Kˣ_J¹::Union{AbstractReciprocalHoppings, Nothing} = nothing,
	Kˣ_J²::Union{AbstractReciprocalHoppings, Nothing} = nothing,
	Kˣ_pair::Union{WanIntijklR, Nothing} = nothing,
	upindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
	dnindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
)
	nw = numorb(Kᵈ_U)
	nw == numorb(Kˣ_U) || error("Mismatched Kᵈ and Kˣ direct interaction!")

	if isnothing(Kᵈ_pair)
		if isnothing(Kᵈ_J¹) && isnothing(Kᵈ_J²)
			W = _U_spinblind_Kᵈ(kgrid, Kᵈ_U)
		else
			if isnothing(Kᵈ_J¹)
				Kᵈ_J¹ = zero(typeof(Kᵈ_J²), nw)
			end
			if isnothing(Kᵈ_J²)
				Kᵈ_J² = zero(typeof(Kᵈ_J¹), nw)
			end
			nw == numorb(Kᵈ_J¹) || error("Mismatched number of wannier basis!")
			nw == numorb(Kᵈ_J²) || error("Mismatched number of wannier basis!")
			W = _U_J_spinblind_Kᵈ(kgrid, Kᵈ_U, Kᵈ_J¹, Kᵈ_J²)
		end
	else
		nw == numorb(Kᵈ_pair) || error("Mismatched number of wannier basis!")
		W = _U_pair_spinblind_Kᵈ(kgrid, Kᵈ_U, Kᵈ_pair)
	end

	if isnothing(Kˣ_pair)
		if isnothing(Kˣ_J¹) && isnothing(Kˣ_J²)
			V = _U_spinblind_Kˣ(kgrid, Kˣ_U)
		else
			if isnothing(Kˣ_J¹)
				Kˣ_J¹ = zero(typeof(Kˣ_J²), nw)
			end
			if isnothing(Kˣ_J²)
				Kˣ_J² = zero(typeof(Kˣ_J¹), nw)
			end
			nw == numorb(Kˣ_J¹) || error("Mismatched number of wannier basis!")
			nw == numorb(Kˣ_J²) || error("Mismatched number of wannier basis!")
			V = _U_J_spinblind_Kˣ(kgrid, Kˣ_U, Kˣ_J¹, Kˣ_J²)
		end
	else
		nw == numorb(Kˣ_pair) || error("Mismatched number of wannier basis!")
		V = _U_pair_spinblind_Kˣ(kgrid, Kˣ_U, Kˣ_pair)
	end

	nk = length(kgrid)
	if isnothing(upindex) && isnothing(dnindex)
		return _spinblind_Kernel(nk, kgrid, W, V)
	else
		if isnothing(upindex) && length(dnindex) == nw
			dnindex = Int.(collect(dnindex))
			upindex = setdiff(1:(nw*2), dnindex)
		elseif isnothing(dnindex) && length(upindex) == nw
			upindex = Int.(collect(upindex))
			dnindex = setdiff(1:(nw*2), upindex)
		else
			error("Wrong spin_index!")
		end
		return _spinaware_Kernel(nk, kgrid, upindex, dnindex, W, V)
	end
end
function Kernel(::Val{:spinaware}, kgrid::RedKgrid,
	Kᵈ_U::DirectInterAction,
	Kˣ_U::DirectInterAction;
	Kᵈ_J¹::Union{AbstractReciprocalHoppings, Nothing} = nothing,
	Kᵈ_J²::Union{AbstractReciprocalHoppings, Nothing} = nothing,
	Kᵈ_pair::Union{WanIntijklR, Nothing} = nothing,
	Kˣ_J¹::Union{AbstractReciprocalHoppings, Nothing} = nothing,
	Kˣ_J²::Union{AbstractReciprocalHoppings, Nothing} = nothing,
	Kˣ_pair::Union{WanIntijklR, Nothing} = nothing,
	upindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
	dnindex::Union{<:AbstractVector{<:Integer}, Nothing} = nothing,
)
	# 只有相互作用包含自旋且需要按照Sz分块计算BSE时才必须使用此情况，其余均使用 `:spinblind` 即可。

	nw = numorb(Kᵈ_U)
	nw == numorb(Kˣ_U) || error("Mismatched Kᵈ and Kˣ interaction!")

	if !isnothing(Kᵈ_J¹) || !isnothing(Kᵈ_J²) || !isnothing(Kˣ_J¹) || !isnothing(Kˣ_J²)
		@info "Don't support J-type interaction with spin. Omit J-type interaction."
	end

	if isnothing(upindex) && isnothing(dnindex)
		error("Have to provide `upindex` or `dnindex`!")
	else
		if isnothing(upindex) && !isnothing(dnindex)
			dnindex = Int.(collect(dnindex))
			upindex = setdiff(1:nw, dnindex)
		elseif !isnothing(upindex) && isnothing(dnindex)
			upindex = Int.(collect(upindex))
			dnindex = setdiff(1:nw, upindex)
		elseif !isnothing(upindex) && !isnothing(dnindex)
			length(union(upindex, dnindex)) == nw || error("Wrong spin_index!")
			upindex = Int.(collect(upindex))
			dnindex = Int.(collect(dnindex))
		else
			error("Wrong spin_index!")
		end
	end

	if isnothing(Kᵈ_pair)
		W = _U_spinaware_Kᵈ(kgrid, Kᵈ_U; upindex, dnindex)
	else
		nw == numorb(Kᵈ_pair) || error("Mismatched number of wannier basis!")
		W = _U_pair_spinaware_Kᵈ(kgrid, Kᵈ_U, Kᵈ_pair; upindex, dnindex)
	end

	if isnothing(Kˣ_pair)
		V = _U_spinaware_Kˣ(kgrid, Kˣ_U; upindex, dnindex)
	else
		nw == numorb(Kˣ_pair) || error("Mismatched number of wannier basis!")
		V = _U_pair_spinaware_Kˣ(kgrid, Kˣ_U, Kˣ_pair; upindex, dnindex)
	end

	nk = length(kgrid)
	return _spinaware_Kernel(nk, kgrid, upindex, dnindex, W, V)
end
