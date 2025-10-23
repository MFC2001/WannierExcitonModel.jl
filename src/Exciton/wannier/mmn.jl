"""
Only need space information, and applicable to any BSE situation.
"""
function _mmn_BSE_uband(bse, uband, λ, bandindex, nnkpts)

	nband = length(bandindex)
	nqb = size(nnkpts, 2)

	tasks = Array{Task}(undef, nband, nband, nqb)
	for i in 1:nqb, β in 1:nband, α in 1:nband
		tasks[α, β, i] = Threads.@spawn begin
			G = ReducedCoordinates(nnkpts[3:5, i])
			if iszero(G)
				M = uband[q].vectors[:, α] ⋅ uband[qb].vectors[:, β]
			else
				G_phase = _uijR_phase(bse, G, λ).phase
				M = uband[q].vectors[:, α] ⋅ (G_phase .* uband[qb].vectors[:, β])
			end
			M
		end
	end
	M = fetch.(tasks)

	return M
end
