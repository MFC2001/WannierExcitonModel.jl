function BSEwannier(qgrid::RedKgrid, bse::BSESU2;
	nnkpfile::Union{AbstractString, Nothing} = nothing,
	outfolder = "./",
)

	if isnothing(nnkpfile)
		nnkpts = kgrid_shell_automatic(qgrid, reciprocal(bse.TB.lattice))
	else
		#TODO
		error("To be continued!")
		nnkpts = Readnnkp(nnkpfile)
	end

	mkpath(outfolder)
	serialize(joinpath(outfolder, "qgrid.sl"), qgrid)
	serialize(joinpath(outfolder, "bse.sl"), bse)
	serialize(joinpath(outfolder, "nnkpts.sl"), nnkpts)

	(BSEband_t, BSEband_u_t, BSEband_s, BSEband_u_s) = BAND(qgrid, bse; vector = true, wfctype = :BlochPeriodic)
	serialize(joinpath(outfolder, "triplet.sl"), BSEband_t)
	serialize(joinpath(outfolder, "singlet.sl"), BSEband_s)
	serialize(joinpath(outfolder, "triplet_u.sl"), BSEband_u_t)
	serialize(joinpath(outfolder, "singlet_u.sl"), BSEband_u_s)

	return nothing
end
function BSEwannier_pp(bsefolder::AbstractString, sym::Symbol;
	outfolder = "./",
	bandindex = nothing,
	guess = nothing,
	guessbasis = :iR,
	λ = 1 // 2,
)

	qgrid = deserialize(joinpath(bsefolder, "qgrid.sl"))
	bse = deserialize(joinpath(bsefolder, "bse.sl"))
	nnkpts = deserialize(joinpath(bsefolder, "nnkpts.sl"))
	if sym == :triplet
		BSEband = deserialize(joinpath(bsefolder, "triplet.sl"))
		BSEband_u = deserialize(joinpath(bsefolder, "triplet_u.sl"))
	elseif sym == :singlet
		BSEband = deserialize(joinpath(bsefolder, "singlet.sl"))
		BSEband_u = deserialize(joinpath(bsefolder, "singlet_u.sl"))
	end

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
	M = _mmn_BSE_uband(bse, BSEband_u, λ, bandindex, nnkpts)
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
