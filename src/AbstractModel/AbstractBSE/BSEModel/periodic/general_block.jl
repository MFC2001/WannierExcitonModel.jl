struct BSEgeneral_block{
	TBT_up <: AbstractTightBindModel,
	TBT_dn <: AbstractTightBindModel,
	KT <: AbstractKernalInterAction,
} <: AbstractBSE
	TB_up::TBT_up
	TB_dn::TBT_dn
	scissor::Float64
	kgrid::RedKgrid
	unitcell::Vector{ReducedCoordinates{Int}}
	vckmap::vckMap
	ijRmap::ijRMap
	bandk::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	Kernal::KT
end
