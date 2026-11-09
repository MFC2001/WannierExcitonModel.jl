
export TightBindModel, numatom, Hamiltonian, spinful

abstract type AbstractTightBindModel <: AbstractModel end

#自旋指标可以由用户控制，wannier基不一定具有确定的自旋。
"""
	TightBindModel{HT<:HermitianReciprocalHoppings} <: AbstractTightBindModel

A tight binding model object containing fields:
- `lattice::Lattice{Float64}`: Lattice of the system.
- `atom_name::Vector{String}`: Names of atoms of the unit cell.
- `atom_location::Vector{ReducedCoordinates{Float64}}`: Locations of atoms in reduced coordinates.
- `orb_name::Vector{String}`: Names of orbitals of the unit cell.
- `orb_location::Vector{ReducedCoordinates{Float64}}`: Locations of orbitals in reduced coordinates.
- `period::Vec3{Bool}`: Periodicity of the system in each lattice basis direction.
- `H::HT`: Hermitian Hamiltonian of the tight binding model.
Note that the auxiliary basis of `lattice` in unperiodic directions should be vertical to the periodic directions.

	numatom(TB::TightBindModel) -> Int
	numorb(TB::TightBindModel) -> Int

Used to obtain the number of atoms or orbitals.
"""
struct TightBindModel{HT<:HermitianReciprocalHoppings} <: AbstractTightBindModel
    lattice::Lattice{Float64}
    atom_name::Vector{String}
    atom_location::Vector{ReducedCoordinates{Float64}}
    orb_name::Vector{String}
    orb_location::Vector{ReducedCoordinates{Float64}}
    period::Vec3{Bool}
    H::HT
end
(TB::TightBindModel)(sym::Symbol, paras...) = TB(Val(sym), paras...)
(TB::TightBindModel)(k::ReducedCoordinates) = TB.H(k)
(TB::TightBindModel)(A::AbstractMatrix, k::ReducedCoordinates) = TB.H(A, k)
(TB::TightBindModel)(k::ReducedCoordinates, orblocat) = TB.H(k, orblocat)
(TB::TightBindModel)(A::AbstractMatrix, k::ReducedCoordinates, orblocat) = TB.H(A, k, orblocat)
(TB::TightBindModel)(::Val{:partial}, k::ReducedCoordinates) =
    TB.H(Val(:partial), k, TB.lattice, TB.orb_location)
(TB::TightBindModel)(::Val{:partial}, A::AbstractArray, k::ReducedCoordinates) =
    TB.H(A::AbstractArray, Val(:partial), k, TB.lattice, TB.orb_location)
function (TB::TightBindModel)(::Val{:orbproj}, aimorbindex::Integer)
    norb = numorb(TB)
    A = zeros(Int, norb, norb)
    A[aimorbindex, aimorbindex] = 1
    return A
end
function (TB::TightBindModel)(::Val{:orbproj}, aimorbindex::AbstractVector{<:Integer})
    norb = numorb(TB)
    A = zeros(Int, norb, norb)
    for i in aimorbindex
        A[i, i] = 1
    end
    return A
end
function (TB::TightBindModel)(::Val{:spinmat}, upindex, dnindex=setdiff(1:numorb(TB), upindex))
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
WannierInterActionBase.numorb(TB::TightBindModel) = length(TB.orb_location)
Hamiltonian(TB::TightBindModel) = TB.H

Base.convert(::TightBindModel{HT₁}, TB::TightBindModel{HT₂}) where {HT₁,HT₂} =
    TightBindModel{HT₁}(TB.lattice, TB.atom_name, TB.atom_location, TB.orb_name, TB.orb_location, TB.period,
        convert(HT₁, TB.H))
Base.convert(::Type{TightBindModel{HT₁}}, TB::TightBindModel{HT₂}) where {HT₁,HT₂} =
    TightBindModel{HT₁}(TB.lattice, TB.atom_name, TB.atom_location, TB.orb_name, TB.orb_location, TB.period,
        convert(HT₁, TB.H))

function Base.show(io::IO, TB::TightBindModel)
    print(io, "$(count(TB.period)) dimensinal Tight binding model with $(numatom(TB)) atoms and $(numorb(TB)) orbitals.")
end

"""
	TightBindModel(; hops::Union{AbstractString,AbstractHR,AbstractReciprocalHoppings},
        cell::Union{AbstractString,Cell},
        orbital::Union{AbstractString,wannier90_centres},
        period=Bool[1, 1, 1],
	)

Create an `TightBindModel` object.
The `AbstractString` is the path to aim file, and check docstrings of `HR`, `AbstractReciprocalHoppings`,
`Cell`, and `wannier90_centres`.
The default format of file is `wannier90_hr`, `POSCAR` and `wannier90_centres`.
"""
function TightBindModel(;
    hops::Union{AbstractString,AbstractHR,AbstractReciprocalHoppings},
    cell::Union{AbstractString,Cell},
    orbital::Union{AbstractString,wannier90_centres},
    period=Bool[1, 1, 1],
)
    if hops isa AbstractString
        hops = ReciprocalHoppings(read(hops, wannier90_hr))
    elseif hops isa AbstractHR
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
function TightBindModel(rh::AbstractReciprocalHoppings, cell::Cell, orbital::wannier90_centres)
    numorb(rh) == numorb(orbital) || error("Mismatched number of orbitals.")
    return TightBindModel(rh, cell, orbital.location; orbname=orbital.name)
end
"""
    TightBindModel(rh::AbstractReciprocalHoppings, cell::Cell, orblocation::AbstractVector{<:AbstractVector{<:Real}};
    orbname::AbstractVector{<:Union{AbstractString,AbstractChar}}=fill("X", length(orblocation)))
"""
function TightBindModel(rh::AbstractReciprocalHoppings, cell::Cell{P₁}, orblocation::AbstractVector{P₂};
    orbname::AbstractVector{<:Union{AbstractString,AbstractChar}}=fill("X", length(orblocation))) where {P₁,P₂<:AbstractVector{<:Real}}
    numorb(rh) == length(orblocation) || error("Mismatched number of orblocation.")
    numorb(rh) == length(orbname) || error("Mismatched number of orbname.")

    if P₁ <: CartesianCoordinates
        cell = convert(cell)
    end
    if P₂ <: ReducedCoordinates
        orblocation = collect(orblocation)
    else
        orblocation = map(CartesianCoordinates{Float64} ∘ collect, orblocation)
        orblocation = collect(cell.lattice \ x for x in orblocation)
    end

    return TightBindModel(
        cell.lattice,
        cell.name,
        cell.location,
        orbname,
        orblocation,
        cell.period,
        HermitianReciprocalHoppings(rh),
    )
end


include("./shell.jl")
# include("./QuantumGeometry/QuantumGeometry.jl")
# include("./BerryCurvature.jl")
# include("./ChernNumber.jl")



function WannierInterActionBase.spinorize(TB::TightBindModel{HT}; mode=conj) where {HT}
    spin_orb_name = [TB.orb_name .* 'u'; TB.orb_name .* 'd']
    spin_orb_location = [TB.orb_location; TB.orb_location]

    spinH = spinorize(TB.H; mode)

    return TightBindModel(
        TB.lattice,
        TB.atom_name,
        TB.atom_location,
        spin_orb_name,
        spin_orb_location,
        TB.period,
        spinH,
    )
end
"""
	merge(TB, hr)
	return a new TB, including hoppings in hr.
"""
function Base.merge(TB::TightBindModel, rh::Union{AbstractHR,AbstractReciprocalHoppings})
    numorb(TB) == numorb(rh) || error("Mismatched orbital.")
    all_rh = merge(TB.H, BaseReciprocalHoppings(rh))
    return TightBindModel(
        deepcopy(TB.lattice),
        copy(TB.atom_name),
        copy(TB.atom_location),
        copy(TB.orb_name),
        copy(TB.orb_location),
        copy(TB.period),
        all_rh,
    )
end
function Base.merge(TB::TightBindModel, arrays::Union{AbstractArray{<:Hopping},Hopping}...)
    all_rh = merge(TB.H, arrays...)
    return TightBindModel(
        deepcopy(TB.lattice),
        copy(TB.atom_name),
        copy(TB.atom_location),
        copy(TB.orb_name),
        copy(TB.orb_location),
        copy(TB.period),
        all_rh,
    )
end
function Base.merge!(TB::TightBindModel, rh::Union{AbstractHR,AbstractReciprocalHoppings})
    numorb(TB) == numorb(rh) || error("Mismatched orbital.")
    merge!(TB.H, BaseReciprocalHoppings(rh))
    return TB
end
function Base.merge!(TB::TightBindModel, arrays::Union{AbstractArray{<:Hopping},Hopping}...)
    merge!(TB.H, arrays...)
    return TB
end
function Base.union(TB::TightBindModel, TBs::TightBindModel...)
    orb_name = vcat(TB.orb_name, [TBi.orb_name for TBi in TBs]...)
    orb_location = vcat(TB.orb_location, [TBi.orb_location for TBi in TBs]...)
    H = union(TB.H, [TBi.H for TBi in TBs]...)
    return TightBindModel(
        TB.lattice,
        TB.atom_name,
        TB.atom_location,
        orb_name,
        orb_location,
        TB.period,
        H,
    )
end
