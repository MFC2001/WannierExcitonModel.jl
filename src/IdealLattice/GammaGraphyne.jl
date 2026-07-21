function γgraphyneTB(; a_intra = 1.425, a_inter = 4.035, c = 20,
	t₀ = 0, t_intra = -2.27, t_inter = 1.8)

	cell = γgraphyneCell(; a_intra, a_inter, c)
	hr = γgraphyneHR(; t₀, t_intra, t_inter)
	orb_location = deepcopy(cell.location)
	orb_name = fill("X", length(orb_location))

	return TightBindModel(
		cell.lattice,
		cell.name,
		cell.location,
		orb_name,
		orb_location,
		HermitianReciprocalHoppings(hr),
		cell.period,
	)
end
