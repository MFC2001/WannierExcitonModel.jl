struct BSEcluster_general_block{
	TBT <: AbstractTightBindModel,
	KT <: KernelInterAction,
} <: AbstractBSE
	TB::TBT
	scissor::Float64
	vcmap::vcMap
	band::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	kernel::KT
end
function Base.show(io::IO, bse::BSEcluster_general_block)
	print(io, "$(count(bse.TB.period)) dimensinal BSE model with $(numatom(bse.TB)) atoms and $(numorb(bse.TB)) orbitals.")
end
