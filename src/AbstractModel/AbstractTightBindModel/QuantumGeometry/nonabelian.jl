"""
	QuantumGeometry(kpoints, rh::AbstractReciprocalHoppings, lattice::Lattice, orblocat::AbstractVector, bandindex)

	Compute the quantum geometry tensor for a set of k-points `kpoints`, given a `ReciprocalHoppings` object `rh`, a `Lattice` object `lattice`, and orbital locations `orblocat`. 
	The optional parameter `bandindex` specifies the band index for which the quantum geometry is computed (default is 1).
"""
function QuantumGeometry(kpoints::AbstractVector{<:ReducedCoordinates}, rh::AbstractReciprocalHoppings, lattice::Lattice, orblocat::AbstractVector, bandindex::AbstractVector{<:Integer})

	norb = numorb(rh)
	nband = length(bandindex)
	nk = length(kgrid)

	QG = Array{ComplexF64}(undef, nband, nband, 3, 3, nk)

	sumindex = setdiff(1:norb, bandindex)

	map(k -> _QG_nonabelian_k!(view(QG, :, :, :, :, k), kpoints[k], rh, lattice, orblocat, sumindex, bandindex), 1:nk)

	return QG
end
function _QG_nonabelian_k!(QG_k, kpoint, rh, lattice, orblocat, sumindex, bandindex)

	nband = length(bandindex)

	E = eigen!(rh(kpoint, orblocat))
	pH = rh(Val(:partial), lattice, kpoint, orblocat)

	for ν in 1:3, μ in 1:3, j in 1:nband, i in 1:nband
		ni = bandindex[i]
		nj = bandindex[j]
		QG_k[i, j, μ, ν] = sum(sumindex) do index
			return (E.vectors[:, ni] ⋅ (pH[:, :, μ] * E.vectors[:, index])) * (E.vectors[:, index] ⋅ (pH[:, :, ν] * E.vectors[:, nj])) /
				   ((E.values[ni] - E.values[index]) * (E.values[nj] - E.values[index]))
		end
	end
	return nothing
end
