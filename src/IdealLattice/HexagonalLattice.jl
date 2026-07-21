function HexagonalTB(; a = 1, c = 20, t₀ = 0, t₁ = 1, t₂ = 0, tₕ = 0)

	cell = HexagonalCell(; a, c)
	hr = HexagonalHR(; t₀, t₁, t₂, tₕ)

	orb_location = deepcopy(cell.location)
	orb_name = fill("X", length(orb_location))

	return TightBindModel(
		cell.lattice,
		cell.name,
		cell.location,
		orb_name,
		orb_location,
		HermitianReciprocalHoppings(hr),
		cell.period)
end
