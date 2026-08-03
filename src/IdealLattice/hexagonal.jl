function IdealTB(::Val{:hexagonal}; a = 1, c = 20, kwargs...)

	cell = IdealCell(Val(:hexagonal); a, c)
	hr = IdealHR(Val(:hexagonal); kwargs...)
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
