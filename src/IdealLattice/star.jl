IdealTB(::Val{:diatomic_kagome}, args...; kwargs...) = IdealTB(Val(:star), args...; kwargs...)
function IdealTB(::Val{:star}; a = 1, b = 2, c = 20, kwargs...)

	cell = IdealCell(Val(:star); a, b, c)
	hr = IdealHR(Val(:star); kwargs...)
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
