
function BSE_kgrid_convergence(TB, v, c, type, nks, nE; q = ReducedCoordinates(0, 0, 0), kwards...)

	eigsolveconfigure!(; howmany = nE, which = EigSorter(abs))

	Gband_t = Matrix{Float64}(undef, nE, length(nks))
	Gband_s = Matrix{Float64}(undef, nE, length(nks))
	for i in eachindex(nks)
		nk = nks[i]
		if iseven(nk)
			kgrid = MonkhorstPack([nk, nk, 1]; kshift = [1 // 2, 1 // 2, 0])
		else
			kgrid = MonkhorstPack([nk, nk, 1])
		end

		bse = BSE(TB, kgrid, v, c, type; kwards...)

		(BSEband_t, BSEband_s) = BAND(q, bse)

		Gband_t[:, i] = BSEband_t
		Gband_s[:, i] = BSEband_s
	end
	return Gband_t, Gband_s

end
