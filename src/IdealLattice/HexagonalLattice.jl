
function HexagonalTB(; a = 1, c = 20, t₀ = 0, t₁ = 1, t₂ = 0, tₕ = 0)

	cell = HexagonalCell(; a, c)
	hr = HexagonalHR(; t₀, t₁, t₂, tₕ)

	orb_location = deepcopy(cell.location)
	orb_name = deepcopy(cell.name)

	return TightBindModel(
		cell.lattice,
		cell.name,
		cell.location,
		orb_name,
		orb_location,
		HermitianReciprocalHoppings(hr),
		cell.period)
end
function HexagonalCell(; a = 1, c = 20)

	l = √3 * a
	𝐚 = [l, 0, 0]
	𝐛 = l * [-1 / 2, √3 / 2, 0]
	𝐜 = [0, 0, c]

	lattice = Lattice([𝐚 𝐛 𝐜])

	location1 = [0, a, c / 2]
	location2 = [l / 2, a / 2, c / 2]

	location = CartesianCoordinates.([location1, location2])
	location = map(x -> lattice \ x, location)

	return Cell(lattice, location; location_type = "ReducedCoordinates", name = fill("X", 2), period = Bool[1, 1, 0])
end
function HexagonalHR(; t₀ = 0, t₁ = 1, t₂ = 0, tₕ = 0)
	allpath = Matrix{Int}(undef, 0, 5)
	allvalue = Vector{Float64}(undef, 0)
	if !iszero(t₀)
		path_t₀ = [
			0 0 0 1 1;
			0 0 0 2 2
		]
		value_t₀ = fill(t₀, size(path_t₀, 1))
		allpath = [allpath; path_t₀]
		allvalue = [allvalue; value_t₀]
	end
	if !iszero(t₁)
		path_t₁ = [
			0 0 0 1 2;
			-1 0 0 1 2;
			0 1 0 1 2;
			0 0 0 2 1;
			1 0 0 2 1;
			0 -1 0 2 1
		]
		value_t₁ = fill(t₁, size(path_t₁, 1))
		allpath = [allpath; path_t₁]
		allvalue = [allvalue; value_t₁]
	end
	if !iszero(t₂)
		path_t₂ = [
			1 1 0 1 1;
			1 0 0 1 1;
			0 1 0 1 1;
			-1 0 0 1 1;
			0 -1 0 1 1;
			-1 -1 0 1 1;
			1 1 0 2 2;
			1 0 0 2 2;
			0 1 0 2 2;
			-1 0 0 2 2;
			0 -1 0 2 2;
			-1 -1 0 2 2
		]
		value_t₂ = fill(t₂, size(path_t₂, 1))
		allpath = [allpath; path_t₂]
		allvalue = [allvalue; value_t₂]
	end
	if !iszero(tₕ)
		path_tₕ = [
			1 1 0 1 1;
			1 0 0 1 1;
			0 1 0 1 1;
			-1 0 0 1 1;
			0 -1 0 1 1;
			-1 -1 0 1 1;
			1 1 0 2 2;
			1 0 0 2 2;
			0 1 0 2 2;
			-1 0 0 2 2;
			0 -1 0 2 2;
			-1 -1 0 2 2
		]
		value_tₕ = tₕ * im * [1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1]
		allpath = [allpath; path_tₕ]
		allvalue = [allvalue; value_tₕ]
	end

	return HR(allpath, allvalue; orbindex = collect(1:2), hrsort = "Y")
end
