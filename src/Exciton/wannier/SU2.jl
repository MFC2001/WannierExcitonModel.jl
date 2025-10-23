function BSEwannier(qgrid::RedKgrid, bse::BSESU2;
	nnkpfile::Union{AbstractString, Nothing} = nothing,
	outfolder = "./",
	bandindex_t = nothing,
	bandindex_s = nothing,
	guess_t = nothing,
	guess_s = nothing,
	guessbasis_t = :iR,
	guessbasis_s = :iR,
	λ = 1 // 2,
	λ_t = λ,
	λ_s = λ,
)

	if isnothing(nnkpfile)
		nnkpts = kgrid_shell_automatic(qgrid, reciprocal(bse.TB.lattice))
	else
		#TODO
		error("To be continued!")
		nnkpts = Readnnkp(nnkpfile)
	end

	bsefolder = joinpath(outfolder, "BSE")
	mkpath(bsefolder)
	serialize(joinpath(bsefolder, "qgrid.sl"), qgrid)
	serialize(joinpath(bsefolder, "bse.sl"), bse)
	serialize(joinpath(bsefolder, "nnkpts.sl"), nnkpts)

	(BSEband_t, BSEband_u_t, BSEband_s, BSEband_u_s) = BAND(qgrid, bse; vector = true, wfctype = :BlochPeriodic)
	serialize(joinpath(bsefolder, "triplet.sl"), BSEband_t)
	serialize(joinpath(bsefolder, "singlet.sl"), BSEband_s)
	serialize(joinpath(bsefolder, "triplet_u.sl"), BSEband_u_t)
	serialize(joinpath(bsefolder, "singlet_u.sl"), BSEband_u_s)

	BSE_wannier_pp(qgrid, bse, nnkpts, BSEband_t, BSEband_u_t;
		outfolder = joinpath(outfolder, "triplet"), bandindex = bandindex_t, guess = guess_t, guessbasis = guessbasis_t, λ = λ_t)
	BSE_wannier_pp(qgrid, bse, nnkpts, BSEband_s, BSEband_u_s;
		outfolder = joinpath(outfolder, "singlet"), bandindex = bandindex_s, guess = guess_s, guessbasis = guessbasis_s, λ = λ_s)

	return nothing
end
function BSE_wannier_pp(qgrid, bse::BSESU2, nnkpts, BSEband, BSEband_u;
	outfolder = "./",
	bandindex = nothing,
	guess = nothing,
	guessbasis = :iR,
	λ = 1 // 2,
)

	allbandindex = collect(1:length(BSEband[1].values))
	if isnothing(bandindex)
		bandindex = allbandindex
	elseif bandindex isa AbstractVector{<:Integer}
		bandindex = Int.(collect(bandindex))
		if bandindex ⊈ allbandindex
			error("bandindex should be a subset of all band index.")
		end
	else
		error("bandindex should be a Vector of Int.")
	end

	BSEband_u = map(uband_getter -> uband_getter(λ), BSEband_u)
	M = _mmn_BSE(qgrid, bse, BSEband_u, bandindex, nnkpts, λ)
	write(joinpath(outfolder, "exciton.mmn"), M, wannier90_mmn; nnkpts, comment = "λ = $λ")
	write(joinpath(outfolder, "exciton.eig"), BSEband, wannier90_eig; bandindex)

	if !isnothing(guess)
		if guessbasis == :iR
			#require guess(ei,Rₑ,hi,Rₕ)
			W = _guess_iR_to_vck(qgrid, bse, guess)
		elseif guessbasis == :ik
			#TODO
		elseif guessbasis == :vck
			#TODO
		else
			error("Wrong wannierbasis.")
		end
		A = _amn_BSE(BSEband, bandindex, W)
		write(joinpath(outfolder, "exciton.amn"), A, wannier90_amn)
	end

	return nothing
end
