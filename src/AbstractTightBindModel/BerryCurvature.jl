function BerryCurvature(::Union{Val{:WilsonLoop}, Val{:WL}},
	kgrid::MonkhorstPack, TB::AbstractTightBindModel; bandindex::Union{Integer, AbstractVector{<:Integer}} = 1, kz::Real = 0)

	BC = BerryCurvature(Val(:WilsonLoop), kgrid, TB.lattice, kz)

	U_vertex = BAND(BC.vertex, TB; vector = true, wfctype = :Periodic)
	U_vertex = map(U -> U.vectors[:, bandindex], U_vertex)
	U_getter = Index_VecOrMat{length(bandindex)}(U_vertex)

	# calculate Berry curvature on each center points
	BC_value = BC(U_getter)

	return BC, BC_value
end
