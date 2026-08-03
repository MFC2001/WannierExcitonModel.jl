struct BSESz{
	TBT_up <: AbstractTightBindModel,
	TBT_dn <: AbstractTightBindModel,
	KT <: KernelInterAction,
} <: AbstractBSE
	TB_up::TBT_up
	TB_dn::TBT_dn
	scissor_c_up::Float64
	scissor_v_up::Float64
	scissor_c_dn::Float64
	scissor_v_dn::Float64
	kgrid::RedKgrid
	unitcell::Vector{ReducedCoordinates{Int}}
	vckmap_uu::vckMap
	vckmap_dd::vckMap
	vckmap_ud::vckMap
	vckmap_du::vckMap
	ijRmap_uu::ijRMap
	ijRmap_dd::ijRMap
	ijRmap_ud::ijRMap
	ijRmap_du::ijRMap
	bandk_up::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	bandk_dn::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	bandkq_up::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	bandkq_dn::Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}
	kernel::KT
	Γdata_up::_BSE_Γdata
	Γdata_dn::_BSE_Γdata
	Huu_upidx::Vector{NTuple{2, Int}}
	Hdd_upidx::Vector{NTuple{2, Int}}
	Huudd_idx::Vector{NTuple{2, Int}}
	Hud_upidx::Vector{NTuple{2, Int}}
	Hdu_upidx::Vector{NTuple{2, Int}}
end
