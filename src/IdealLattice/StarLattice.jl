function StarTB(; a = 1, b = 2, c = 20,
	t₀ = 0, t₁ = 1.5, t₂ = 1, t₃ = 0, t₄ = 0, t₅ = 0, t_inverse = 0)

	cell = StarCell(; a, b, c)
	orb_location = deepcopy(cell.location)
	orb_name = fill("X", length(orb_location))

	hr = StarHR(; t₀, t₁, t₂, t₃, t₄, t₅, t_inverse)

	return TightBindModel(
		cell.lattice,
		cell.name,
		cell.location,
		orb_name,
		orb_location,
		HermitianReciprocalHoppings(hr),
		cell.period)
end
