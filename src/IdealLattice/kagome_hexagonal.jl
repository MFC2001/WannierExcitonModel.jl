function IdealTB(::Val{:kagome_hexagonal}; a = 1, b = 2, c = 20, kwargs...)

	cell = IdealCell(Val(:kagome_hexagonal); a, b, c)
	hr = IdealHR(Val(:kagome_hexagonal); kwargs...)
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
