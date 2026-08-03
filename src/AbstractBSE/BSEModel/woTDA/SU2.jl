struct BSESU2_woTDA{
	eleHT <: HermitianReciprocalHoppings,
	KT <: KernelInterAction,
} <: AbstractBSE
	lattice::Lattice{Float64}
	atom_name::Vector{String}
	atom_location::Vector{ReducedCoordinates{Float64}}
	orb_name::Vector{String}
	orb_location::Vector{ReducedCoordinates{Float64}}
	period::Vec3{Bool}
	eleH::eleHT
	scissor::Float64
	kgrid::RedKgrid
	unitcell::Vector{ReducedCoordinates{Int}}
	vckmap::vckMap
	ijRmap::ijRMap
	bandk::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	bandkq::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	kernel::KT
	Γdata::_BSE_Γdata
	H_upidx::Vector{NTuple{2, Int}}
end
