function BSE_wannier(qgrid::RedKgrid, bse::BSEspinless;
	nnkpfile::Union{AbstractString, Nothing} = nothing,
	outfolder = "./",
	bandindex_t = nothing,
	bandindex_s = nothing,
	guess_t = nothing,
	guess_s = nothing,
	guessbasis_t = :iR,
	guessbasis_s = :iR,
	η = 1 // 2,
	η_t = η,
	η_s = η,
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

	(BSEband_t, BSEband_u_t, BSEband_s, BSEband_u_s) = BAND(qgrid, bse; vector = true, wfctype = :BlochPeriodic, η_t, η_s)
	serialize(joinpath(bsefolder, "triplet.sl"), BSEband_t)
	serialize(joinpath(bsefolder, "singlet.sl"), BSEband_s)

	BSE_wannier_pp(qgrid, bse, nnkpts, BSEband_t;
		outfolder = joinpath(outfolder, "triplet"), bandindex = bandindex_t, guess = guess_t, guessbasis = guessbasis_t, η = η_t)
	BSE_wannier_pp(qgrid, bse, nnkpts, BSEband_s;
		outfolder = joinpath(outfolder, "singlet"), bandindex = bandindex_s, guess = guess_s, guessbasis = guessbasis_s, η = η_s)

	return nothing
end
function BSE_wannier_pp(qgrid, bse::BSEspinless, nnkpts, BSEband;
	outfolder = "./",
	bandindex = nothing,
	guess = nothing,
	guessbasis = :iR,
	η = 1 // 2,
)

	if isnothing(bandindex)
		bandindex = collect(1:length(BSEband[1].values))
	elseif bandindex isa AbstractVector{<:Integer}
		bandindex = Int.(collect(bandindex))
	else
		error("bandindex should be a Vector of Int.")
	end

	M = _mmn_BSE(qgrid, bse, BSEband, bandindex, nnkpts, η)
	Writemmn(nnkpts, M, joinpath(outfolder, "exciton.mmn"))
	Writeeig(BSEband, joinpath(outfolder, "exciton.eig"); bandindex)

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
		Writeamn(A, joinpath(outfolder, "exciton.amn"))
	end

	return nothing
end
