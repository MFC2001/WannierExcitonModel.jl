function IdealTB(::Val{:γgraphyne}; a_intra = 1.425, a_inter = 4.035, c = 20, kwargs...)

	cell = IdealCell(Val(:γgraphyne); a_intra, a_inter, c)
	hr = IdealHR(Val(:γgraphyne); kwargs...)
	orb_location = deepcopy(cell.location)
	orb_name = fill("X", length(orb_location))

	return TightBindModel(
		cell.lattice,
		cell.name,
		cell.location,
		orb_name,
		orb_location,
		cell.period,
		HermitianReciprocalHoppings(hr),
	)
end
