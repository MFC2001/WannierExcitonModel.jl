
export TightBindModel

abstract type AbstractTightBindModel <: AbstractModel end

#自旋指标可以由用户控制，wannier基不一定具有确定的自旋。
"""
	TightBindModel{T <: Number, U <: AbstractReciprocalHoppings{T}} <: AbstractTightBindModel

A tight binding model object containing fields:
- `lattice::Lattice{Float64}`: Lattice of the system.
- `atom_name::Vector{String}`: Names of atoms of the unit cell.
- `atom_location::Vector{ReducedCoordinates{Float64}}`: Locations of atoms in reduced coordinates.
- `orb_name::Vector{String}`: Names of orbitals of the unit cell.
- `orb_location::Vector{ReducedCoordinates{Float64}}`: Locations of orbitals in reduced coordinates.
- `H::HermitianReciprocalHoppings{T, U}`: Hermitian Hamiltonian of the tight binding model.
- `period::Vec3{Bool}`: Periodicity of the system in each lattice basis direction.
Note that the auxiliary basis of `lattice` in unperiodic directions should be vertical to the periodic directions.

	numatom(TB::TightBindModel) -> Int
	numorb(TB::TightBindModel) -> Int

Used to obtain the number of atoms or orbitals.
"""
struct TightBindModel{T <: Number, U <: AbstractReciprocalHoppings{T}} <: AbstractTightBindModel
	lattice::Lattice{Float64}
	atom_name::Vector{String}
	atom_location::Vector{ReducedCoordinates{Float64}}
	orb_name::Vector{String}
	orb_location::Vector{ReducedCoordinates{Float64}}
	H::HermitianReciprocalHoppings{T, U}
	period::Vec3{Bool}
end
(TB::TightBindModel)(sym::Symbol, paras...) = TB(Val(sym), paras...)
(TB::TightBindModel)(k::ReducedCoordinates) = TB.H(k)
(TB::TightBindModel)(A::AbstractMatrix, k::ReducedCoordinates) = TB.H(A, k)
(TB::TightBindModel)(k::ReducedCoordinates, orblocat) = TB.H(k, orblocat)
(TB::TightBindModel)(A::AbstractMatrix, k::ReducedCoordinates, orblocat) = TB.H(A, k, orblocat)
(TB::TightBindModel)(::Val{:partial}, k::ReducedCoordinates) =
	TB.H(Val(:partial), TB.lattice, k, TB.orb_location)
(TB::TightBindModel)(::Val{:partial}, A::AbstractArray, k::ReducedCoordinates) =
	TB.H(A::AbstractArray, Val(:partial), TB.lattice, k, TB.orb_location)

function (TB::TightBindModel)(::Val{:spinmat}, upindex, dnindex = setdiff(1:numorb(TB), upindex))
	norb = numorb(TB)
	A = zeros(Int, norb, norb)
	for i in upindex
		A[i, i] = 1
	end
	for i in dnindex
		A[i, i] = -1
	end
	return A
end

numatom(TB::TightBindModel) = length(TB.atom_location)
numorb(TB::TightBindModel) = length(TB.orb_location)
Hamiltonian(TB::TightBindModel) = TB.H

Base.convert(::TightBindModel{T₁, U₁}, TB::TightBindModel{T₂, U₂}) where {T₁, U₁, T₂, U₂} =
	TightBindModel{T₁, U₁}(TB.lattice, TB.atom_name, TB.atom_location, TB.orb_name, TB.orb_location,
		convert(HermitianReciprocalHoppings{T₁, U₁}, TB.H), TB.period)
Base.convert(::Type{TightBindModel{T₁, U₁}}, TB::TightBindModel{T₂, U₂}) where {T₁, U₁, T₂, U₂} =
	TightBindModel{T₁, U₁}(TB.lattice, TB.atom_name, TB.atom_location, TB.orb_name, TB.orb_location,
		convert(HermitianReciprocalHoppings{T₁, U₁}, TB.H), TB.period)

function Base.show(io::IO, TB::TightBindModel)
	print(io, "$(count(TB.period)) dimensinal Tight binding model with $(numatom(TB)) atoms and $(numorb(TB)) orbitals.")
end

"""
	TightBindModel(;
		hops::Union{AbstractString, HR, AbstractReciprocalHoppings},
		cell::Union{AbstractString, Cell},
		orbital::Union{AbstractString, ORBITAL},
		period = Bool[1, 1, 1],
	)

Create an `TightBindModel` object.
The `AbstractString` is the path to aim file, and check docstrings of `HR`, `AbstractReciprocalHoppings`, `Cell` and `ORBITAL`.
"""
function TightBindModel(;
	hops::Union{AbstractString, HR, AbstractReciprocalHoppings},
	cell::Union{AbstractString, Cell},
	orbital::Union{AbstractString, ORBITAL},
	period = Bool[1, 1, 1],
)
	if hops isa AbstractString
		hops = ReciprocalHoppings(read(hops, wannier90_hr))
	elseif hops isa HR
		hops = ReciprocalHoppings(hops)
	end
	if cell isa AbstractString
		cell = read(cell, POSCAR; period)
	end
	if orbital isa AbstractString
		orbital = read(orbital, wannier90_centres)
	end
	return TightBindModel(hops, cell, orbital)
end
function TightBindModel(rh::AbstractReciprocalHoppings, cell::Cell{<:CartesianCoordinates}, orbital::ORBITAL)
	return TightBindModel(rh, convert(cell), orbital)
end
function TightBindModel(rh::AbstractReciprocalHoppings, cell::Cell{<:ReducedCoordinates}, orbital::ORBITAL)
	numorb(rh) == numorb(orbital) || error("Mismatched number of orbitals.")
	return TightBindModel(
		cell.lattice,
		cell.name,
		cell.location,
		orbital.name,
		[cell.lattice \ x for x in orbital.location],
		HermitianReciprocalHoppings(rh),
		cell.period,
	)
end

include("./shell.jl")
# include("./QuantumGeometry/QuantumGeometry.jl")
# include("./BerryCurvature.jl")
# include("./ChernNumber.jl")





function spinTightBindModel(TB::TightBindModel{T, U}; mode = conj) where {T, U}
	spin_orb_name = [TB.orb_name .* 'u'; TB.orb_name .* 'd']
	spin_orb_location = [TB.orb_location; TB.orb_location]

	spinH = spinrh(TB.H; mode = conj)

	return TightBindModel(
		TB.lattice,
		TB.atom_name,
		TB.atom_location,
		spin_orb_name,
		spin_orb_location,
		spinH,
		TB.period,
	)
end
"""
	union(TB, hr)
	return a new TB, including hoppings in hr.
"""
function Base.union(TB::TightBindModel{T₁, U}, hr::HR{T₂}) where {T₁, U, T₂}
	numorb(TB) == numorb(hr) || error("Mismatched orbital.")
	rh = ReciprocalHoppings(hr)
	return union(TB, rh)
end
function Base.union(TB::TightBindModel{T₁, U}, rh::AbstractReciprocalHoppings{T₂}) where {T₁, U, T₂}
	all_rh = union(TB.H, rh)
	return TightBindModel(
		deepcopy(TB.lattice),
		copy(TB.atom_name),
		copy(TB.atom_location),
		copy(TB.orb_name),
		copy(TB.orb_location),
		HermitianReciprocalHoppings(all_rh),
		copy(TB.period),
	)
end



