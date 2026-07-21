function absorption(bse::Union{BSESU2, BSEgeneral}, ωgrid::AbstractVector{<:Real};
	ê = [0, 0, 1], broadening::Symbol = :gauss, scale::Real = 0.1, outfolder = nothing,
)

	nk = length(bse.kgrid)
	Ω = abs(det(parent(bse.TB.lattice)))
	coff = π * (qₑ * 1000.0) / (ϵ₀ * nk * Ω)

	Nvck = length(bse.vckmap)

	ΔE = Vector{Float64}(undef, Nvck)
	eupHu2 = Vector{Float64}(undef, Nvck) # v = upHu / ħ
	Threads.@threads for i in Base.OneTo(Nvck)
		(v, c, k) = bse.vckmap[i]
		ΔE[i] = bse.bandk[k].values[c] - bse.bandk[k].values[v] # eV
		eupHu2[i] = abs2(ê ⋅ bse.Γdata.upHu[:, i]) / ΔE[i]^2 # Å^2
	end

	result = broaden_peaks(broadening, ωgrid, ΔE, eupHu2; scale, normalize = true)
	result .*= coff

	if !isnothing(outfolder)
		mkpath(outfolder)
		_write_absorption(joinpath(outfolder, "absorption_noeh.dat"), ωgrid, result)
	end
	return result
end
function absorption(bse::Union{BSESU2, BSEgeneral}, ωgrid::AbstractVector{<:Real}, eigen_vck::Eigen;
	ê = [0, 0, 1], broadening::Symbol = :gauss, scale::Real = 0.1, outfolder = nothing,
)

	nk = length(bse.kgrid)
	Ω = abs(det(parent(bse.TB.lattice)))
	coff = π * (qₑ * 1000.0) / (ϵ₀ * nk * Ω)

	Nvck = length(bse.vckmap)

	ΔE = Vector{Float64}(undef, Nvck)
	eupHu = Vector{ComplexF64}(undef, Nvck) # eV ⋅ Å, v = upHu / ħ
	Threads.@threads for i in Base.OneTo(Nvck)
		(v, c, k) = bse.vckmap[i]
		ΔE[i] = bandk[k].values[c] - bandk[k].values[v] # eV
		eupHu[i] = ê ⋅ bse.Γdata.upHu[:, i]
	end
	eupHu2 = map(i->abs2(eupHu[i]) / ΔE[i]^2, Base.OneTo(Nvck)) # Å^2

	Nexc = length(eigen_vck.values)
	AeupHu2 = Vector{Float64}(undef, Nexc) # Å^2
	for i in 1:Nexc
		AeupHu2[i] = abs2(sum(ivck->eigen_vck.vectors[ivck, i] * eupHu[ivck], Base.OneTo(Nvck))) / eigen_vck.values[i]^2
	end

	result_woeh = broaden_peaks(broadening, ωgrid, ΔE, eupHu2; scale, normalize = true)
	result_wieh = broaden_peaks(broadening, ωgrid, eigen_vck.values, AeupHu2; scale, normalize = true)
	result_woeh .*= coff
	result_wieh .*= coff

	if !isnothing(outfolder)
		mkpath(outfolder)
		_write_absorption(joinpath(outfolder, "absorption_noeh.dat"), ωgrid, result_woeh)
		_write_absorption(joinpath(outfolder, "absorption_eh.dat"), ωgrid, result_wieh)
	end
	return result_wieh, result_woeh
end
