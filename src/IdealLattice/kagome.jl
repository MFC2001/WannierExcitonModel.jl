function IdealTB(::Val{:kagome}; a = 1, b = a, c = 20, kwargs...)

	cell = IdealCell(Val(:kagome); a, b, c)
	hr = IdealHR(Val(:kagome); kwargs...)
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
