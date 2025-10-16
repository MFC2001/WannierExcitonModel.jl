export fold2FBZ

const _neighbour_Γ = [
	[-1, -1, -1],
	[-1, -1, 0],
	[-1, -1, 1],
	[-1, 0, -1],
	[-1, 0, 0],
	[-1, 0, 1],
	[-1, 1, -1],
	[-1, 1, 0],
	[-1, 1, 1],
	[0, -1, -1],
	[0, -1, 0],
	[0, -1, 1],
	[0, 0, -1],
	[0, 0, 1],
	[0, 1, -1],
	[0, 1, 0],
	[0, 1, 1],
	[1, -1, -1],
	[1, -1, 0],
	[1, -1, 1],
	[1, 0, -1],
	[1, 0, 0],
	[1, 0, 1],
	[1, 1, -1],
	[1, 1, 0],
	[1, 1, 1],
]

"""
	fold2FBZ(Rlattice::ReciprocalLattice, kpoints::AbstractVector{<:AbstractVector}) -> AbstractVector{<:AbstractVector}
	fold2FBZ(Rlattice::ReciprocalLattice, kpoint::AbstractVector{<:Real}) -> AbstractVector{<:Real}

Return the input's equivalent kpoints in FBZ.
"""
function fold2FBZ(Rlattice::ReciprocalLattice, kpoints::AbstractVector{<:AbstractVector})

	_neighbour_Γ_car = map(nΓ -> Rlattice * nΓ, _neighbour_Γ)

	foldedkpoints = deepcopy(kpoints)
	for kindex in eachindex(foldedkpoints)

		kpoint = foldedkpoints[kindex]

		if _is_approx_integer(kpoint; atol = 1e-4)
			foldedkpoints[kindex] = kpoint - round.(kpoint)
			continue
		end

		kpoint_coord = Rlattice * kpoint

		notinFBZ = true
		while notinFBZ
			Δ2Γ = norm(kpoint_coord)
			notinFBZ = false
			for (i, nΓ) in enumerate(_neighbour_Γ)
				T = kpoint_coord - _neighbour_Γ_car[i]
				if norm(T) < Δ2Γ
					kpoint_coord = T
					kpoint = kpoint - nΓ
					notinFBZ = true
					break
				end
			end
		end

		foldedkpoints[kindex] = kpoint
	end

	return foldedkpoints
end

function fold2FBZ(Rlattice::ReciprocalLattice, kpoint::AbstractVector{<:Real})

	if _is_approx_integer(kpoint; atol = 1e-6)
		return kpoint - round.(kpoint)
	end

	kpoint_coord = Rlattice * kpoint
	_neighbour_Γ_car = map(nΓ -> Rlattice * nΓ, _neighbour_Γ)

	notinFBZ = true
	while notinFBZ
		Δ2Γ = norm(kpoint_coord)
		notinFBZ = false
		for (i, nΓ) in enumerate(_neighbour_Γ)
			T = kpoint_coord - _neighbour_Γ_car[i]
			if norm(T) < Δ2Γ
				kpoint_coord = T
				kpoint = kpoint - nΓ
				notinFBZ = true
				break
			end
		end
	end

	return kpoint
end
