
function BerryCurvature(::Union{Val{:WilsonLoop}, Val{:WL}},
	kgrid::MonkhorstPack, bse::BSEspinless; bandindex::Union{Integer, AbstractVector{<:Integer}} = 1, kz::Real = 0)

	BC = BerryCurvature(Val(:WilsonLoop), kgrid, bse.TB.lattice, kz)

	bandt_vertex, bands = BAND(BC.vertex, bse; vector = true, wfctype = :Periodic)
	U_vertex = map(U -> U.vectors[:, bandindex], U_vertex)
	U_getter = _BC_ugetter{length(bandindex)}(U_vertex)

	# calculate Berry curvature on each center points
	BC_value = BC(U_getter)

	return BC, BC_value
end



# function BerryCurvature_pre(::Union{Val{:WilsonLoop}, Val{:WL}}, kgrid::MonkhorstPack, bse::BSEspinless, kz::Real = 0)

# 	center, vertex, center_vertex_loop = _BerryCurvature_center_vertex(kgrid.kgrid_size, kz)
# 	ΔS = _BerryCurvature_microsurface(kgrid.kgrid_size, TB.lattice)

# 	# precalculation of U matrices on vertices
# 	U_vertex = BAND(vertex, bse; vector = true)

# 	# fold all kpoints into the first Brillouin zone
# 	center = fold2FBZ(reciprocal(TB.lattice), center)
# 	kshift = map(nk -> isodd(nk) ? 0 : 1 // 2, kgrid.kgrid_size)
# 	center = RedKgrid(center, kgrid.kgrid_size, kshift)

# 	return center, center_vertex_loop, U_vertex, ΔS
# end
# function BerryCurvature(::Union{Val{:WilsonLoop}, Val{:WL}},
# 	kgrid::MonkhorstPack, bse::BSEspinless, bandindex::Union{Integer, AbstractVector{<:Integer}}, kz::Real = 0)

# 	center, center_vertex_loop, U_vertex, ΔS = BerryCurvature_pre(Val(:WilsonLoop), kgrid, TB, kz)

# 	U_vertex = map(U -> U.vectors[:, bandindex], U_vertex)
# 	U_getter = _BC_ugetter{length(bandindex)}(U_vertex)

# 	# calculate Berry curvature on each center points
# 	BC = BerryCurvature(Val(:WilsonLoop), U_getter, center_vertex_loop, ΔS)

# 	return center, BC
# end
