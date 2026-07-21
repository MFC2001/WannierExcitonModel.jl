
function QuantumGeometry(kpoints::AbstractVector{<:ReducedCoordinates}, rh::AbstractReciprocalHoppings, lattice::Lattice, orblocat::AbstractVector, bandindex::Integer)

	nk = length(kpoints)

	QG = Array{ComplexF64}(undef, 3, 3, nk)

	sumindex = setdiff(1:numorb(rh), bandindex)

	map(k -> _QG_abelian_k!(view(QG, :, :, k), kpoints[k], rh, lattice, orblocat, sumindex, bandindex), 1:nk)

	return QG
end
function _QG_abelian_k!(QG_k, kpoint, rh, lattice, orblocat, sumindex, bandindex)

	E = eigen!(rh(kpoint, orblocat))
	pH = rh(Val(:partial), lattice, kpoint, orblocat)

	E₀ = E.values[bandindex]
	φ₀ = E.vectors[:, bandindex]
	for μ in 1:3, ν in 1:3
		QG_k[μ, ν] = sum(sumindex) do index
			return (φ₀ ⋅ (pH[:, :, μ] * E.vectors[:, index])) * (E.vectors[:, index] ⋅ (pH[:, :, ν] * φ₀)) / (E₀ - E.values[index])^2
		end
	end

	return nothing
end
