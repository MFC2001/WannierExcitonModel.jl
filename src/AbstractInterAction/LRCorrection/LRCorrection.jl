
export UwithLR, AbstractLRCorrection

"""
	abstract type AbstractLRCorrection end

Long-range correction for coulomb interaction, can give the potential matrix in reciprocal space.
It supports 

	(v::AbstractLRCorrection)(k::AbstractVector)
	(v::AbstractLRCorrection)(k::AbstractVector, nk::Integer)
	(v::AbstractLRCorrection)(A::AbstractMatrix, k::AbstractVector)
	(v::AbstractLRCorrection)(A::AbstractMatrix, k::AbstractVector, nk::Integer)
	(v::AbstractLRCorrection)(::Val{+}, A, k::AbstractVector)
	(v::AbstractLRCorrection)(::Val{+}, A, k::AbstractVector, nk::Integer)

The methods involving `nk` can give head term when k is zero in all periodic directions.
"""
abstract type AbstractLRCorrection end

include("./GaussLRCorrection.jl")

# We calculate head term instead of k=0.
# 这里的0D长程修正仅用于不考虑长程修正的情况。
# 默认在有长程修正的情况下输出结果一定是复的。
"""
	UwithLR{<: AbstractReciprocalHoppings, <: AbstractLRCorrection} <: AbstractInterAction

Direct coulomb term between wannier basis with long-range correction.
Fields:
- `norb::Int`: the number of orbitals;
- `SR::AbstractReciprocalHoppings`: short-range part of direct term;
- `LR::AbstractLRCorrection`: long-range part of direct term.

We have achieved long-range correction by using Gauss potential.
If you want to use other potential, you can define a new subtype of `AbstractLRCorrection`. 
Make sure your subtype has methods:

	(v::yoursubtype)(::Val{+}, A, k::AbstractVector)
	(v::yoursubtype)(::Val{+}, A, k::AbstractVector, nk::Integer)

For an instance `u` of `UwithLR`, you can run:

```julia
julia> u(k)
julia> u(k, nk)
```

These two methods will create a new matrix. And `u(k, nk)` will judge whether head term or normal value by k.

```julia
julia> u(A, k)
julia> u(A, k, nk)
```

These two methods will change `A` instead of creating a new matrix.
"""
struct UwithLR{S <: AbstractReciprocalHoppings, L <: AbstractLRCorrection} <: AbstractInterAction
	norb::Int
	SR::S
	LR::L
end
function (v::UwithLR)(k::AbstractVector)
	A = Matrix{ComplexF64}(undef, v.norb, v.norb)
	return v(A, k)
end
function (v::UwithLR)(k::AbstractVector, nk::Integer)
	A = Matrix{ComplexF64}(undef, v.norb, v.norb)
	return v(A, k, nk)
end
function (v::UwithLR)(A::AbstractMatrix, k::AbstractVector)
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	v.SR(A, k)
	v.LR(Val(+), A, k)
	return A
end
function (v::UwithLR)(A::AbstractMatrix, k::AbstractVector, nk::Integer)
	size(A) == (v.norb, v.norb) || error("Buffer size mismatch.")
	v.SR(A, k)
	v.LR(Val(+), A, k, nk)
	return A
end

function Base.show(io::IO, U::UwithLR{S, L}) where {S, L <: AbstractLRCorrection}
	println(io, "Direct term with long-range correction between electrionic wannier bases.")
end

"""
	UwithLR(U, lattice::Lattice, orblocat::AbstractVector{<:ReducedCoordinates}, rcut::Real; kwargs...) -> UwithLR

Construct a `UwithLR` instance with given short-range part `U`, lattice, orbital locations and cutoff radius `rcut`.
This long-range correction is achieved by using Gauss potential.

- `U`: the direct Coulomb potential between the wannier bases within a limited distance, can be a `path::String`, `hr::HR` or `rh::AbstractReciprocalHoppings`.
kwargs:
- `period`: the periodicity;
- `ϵ::Real`: dielectric constant;
- `αrcut::Real`: α × rcut, default is 4.5;
- `δ:Real`: the cutoff value of reciprocal potential.

Make sure the input U has converged to 1/r.
`rcut` should be less than a real radius value, out of which the U approximate to 1/r.

!!! note
	The less rcut is, the more kpoints needed to calculate long-range correction, so use a rcut as larger as possible.
"""
function UwithLR(U::AbstractString, paras...; kwards...)
	U = ReciprocalHoppings(read(U, wannier90_hr))
	return UwithLR(U, paras...; kwards...)
end
function UwithLR(U::HR, paras...; kwards...)
	U = ReciprocalHoppings(U)
	return UwithLR(U, paras...; kwards...)
end
function UwithLR(U::AbstractReciprocalHoppings, lattice::Lattice, orblocat::AbstractVector{<:ReducedCoordinates}, rcut::Real;
	period::AbstractVector{<:Union{Bool, Integer}} = Bool[1, 1, 1], αrcut::Real = 4.5, δ::Real = 1e-6, ϵ::Real = 1)

	numorb(U) == length(orblocat) || error("Wrong number of orbital locations.")
	length(period) == 3 || error("Wrong period.")

	period = Vec3{Bool}(collect(period))
	p = count(period)
	if p == 0
		#建议为这种情况使用另一种方案，不将其定义为UwithLR.
		return _UwithLR_0D(U, lattice, orblocat, ϵ)
	end

	α = αrcut / rcut

	V_SR = _UwithLR_SR(lattice, orblocat, U, α, ϵ)
	if p == 3
		V_LR = GaussLRCorrection3D(lattice, orblocat, α; δ, ϵ)
	elseif p == 2
		V_LR = GaussLRCorrection2D(lattice, orblocat, period, α; δ, ϵ)
	elseif p == 1
		#TODO 
		error("TODO")
		# V_LR = GaussLRCorrection1D(lattice, orblocat, period, α; δ, ϵ)
	end

	return UwithLR(length(orblocat), V_SR, V_LR)
end
function UwithLR(SR::AbstractReciprocalHoppings{T}) where {T}
	#This means you don't want include long-range correction.
	SR = ReciprocalHoppings(SR)
	LR = GaussLRCorrection0D(numorb(SR))
	return UwithLR{T, GaussLRCorrection0D}(numorb(SR), SR, LR)
end

function _UwithLR_SR(lattice, orblocat, U, α, ϵ)

	φR = RealGauss(; ϵ, α)

	SR = deepcopy(U)
	hops_SR = hops(SR)

	for I in CartesianIndices(hops_SR)
		(i, j) = Tuple(I)
		for (ii, hop) in enumerate(hops_SR[I])
			value_SR = hop.t - φR(lattice * (hop.R + orblocat[j] - orblocat[i]))
			hops_SR[I][ii] = similar(hop, value_SR)
		end
	end

	return SR
end

function _UwithLR_0D(U::AbstractReciprocalHoppings{UT}, lattice, orblocat, ϵ) where {UT}

	SR = deepcopy(U)
	hops_SR = hops(SR)
	Nhop_SR = Nhop(SR)

	I = findfirst(n -> n > 1, Nhop_SR)
	if !isnothing(I)
		error("Wrong 0D U.")
	end

	VR = RealInverseR(; ϵ)

	I = findall(iszero, Nhop_SR)
	for II in I
		(i, j) = Tuple(II)
		push!(hops_SR[II], Hopping{UT}(i, j, [0, 0, 0], VR(lattice * (orblocat[j] - orblocat[i]))))
		Nhop_SR[II] = 1
	end

	norb = numorb(SR)
	LR = GaussLRCorrection0D(norb)

	return UwithLR(norb, SR, LR)
	# return SR
end

