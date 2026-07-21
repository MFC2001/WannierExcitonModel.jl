function absorption(bse::BSESz, ωgrid::AbstractVector{<:Real};
	ê = [0, 0, 1], broadening::Symbol = :gauss, scale::Real = 0.1, outfolder = nothing,
)

	nk = length(bse.kgrid)
	Ω = abs(det(parent(bse.TB.lattice)))
	coff = π * (qₑ * 1000.0) / (ϵ₀ * nk * Ω)

	Nvck_uu = length(bse.vckmap_uu)
	if Nvck_uu == 0
		result_uu = zeros(Float64, length(ωgrid))
	else
		ΔE = Vector{Float64}(undef, Nvck_uu)
		eupHu2 = Vector{Float64}(undef, Nvck_uu) # v = upHu / ħ
		Threads.@threads for i in Base.OneTo(Nvck_uu)
			(v, c, k) = bse.vckmap_uu[i]
			ΔE[i] = bse.bandk_up[k].values[c] - bse.bandk_up[k].values[v] # eV
			eupHu2[i] = abs2(ê ⋅ bse.Γdata_up.upHu[:, i]) / ΔE[i]^2 # Å^2
		end

		result_uu = broaden_peaks(broadening, ωgrid, ΔE, eupHu2; scale, normalize = true)
		result_uu .*= coff
	end

	Nvck_dd = length(bse.vckmap_dd)
	if Nvck_dd == 0
		result_dd = zeros(Float64, length(ωgrid))
	else
		ΔE = Vector{Float64}(undef, Nvck_dd)
		eupHu2 = Vector{Float64}(undef, Nvck_dd) # v = upHu / ħ
		Threads.@threads for i in Base.OneTo(Nvck_dd)
			(v, c, k) = bse.vckmap_dd[i]
			ΔE[i] = bse.bandk_dn[k].values[c] - bse.bandk_dn[k].values[v] # eV
			eupHu2[i] = abs2(ê ⋅ bse.Γdata_dn.upHu[:, i]) / ΔE[i]^2 # Å^2
		end

		result_dd = broaden_peaks(broadening, ωgrid, ΔE, eupHu2; scale, normalize = true)
		result_dd .*= coff
	end

	if !isnothing(outfolder)
		mkpath(outfolder)
		_write_absorption(joinpath(outfolder, "absorption_noeh_up.dat"), ωgrid, result_uu)
		_write_absorption(joinpath(outfolder, "absorption_noeh_dn.dat"), ωgrid, result_dd)
	end
	return result_uu, result_dd
end
function absorption(bse::BSESz, ωgrid::AbstractVector{<:Real}, eigen_vck::Eigen;
	ê = [0, 0, 1], broadening::Symbol = :gauss, scale::Real = 0.1, outfolder = nothing,
)

	# eigen_vck should have Sz = 0

	nk = length(bse.kgrid)
	Ω = abs(det(parent(bse.TB.lattice)))
	coff = π * (qₑ * 1000.0) / (ϵ₀ * nk * Ω)

	Nexc = length(eigen_vck.values)
	AeupHu = zero(ComplexF64, Nexc) # Å^2

	Nvck_uu = length(bse.vckmap_uu)
	if Nvck_uu == 0
		result_uu = zeros(Float64, length(ωgrid))
	else
		ΔE = Vector{Float64}(undef, Nvck_uu)
		eupHu = Vector{Float64}(undef, Nvck_uu) # eV ⋅ Å, v = upHu / ħ
		Threads.@threads for i in Base.OneTo(Nvck_uu)
			(v, c, k) = bse.vckmap_uu[i]
			ΔE[i] = bse.bandk_up[k].values[c] - bse.bandk_up[k].values[v] # eV
			eupHu[i] = ê ⋅ bse.Γdata_up.upHu[:, i]
		end
		eupHu2 = map(i->abs2(eupHu[i]) / ΔE[i]^2, Base.OneTo(Nvck_uu)) # Å^2
		result_uu = broaden_peaks(broadening, ωgrid, ΔE, eupHu2; scale, normalize = true)
		result_uu .*= coff

		for i in 1:Nexc
			AeupHu[i] += sum(ivck->eigen_vck.vectors[ivck, i] * eupHu[ivck], Base.OneTo(Nvck_uu))
		end
	end

	Nvck_dd = length(bse.vckmap_dd)
	if Nvck_dd == 0
		result_dd = zeros(Float64, length(ωgrid))
	else
		ΔE = Vector{Float64}(undef, Nvck_dd)
		eupHu = Vector{Float64}(undef, Nvck_dd) # eV ⋅ Å, v = upHu / ħ
		Threads.@threads for i in Base.OneTo(Nvck_dd)
			(v, c, k) = bse.vckmap_dd[i]
			ΔE[i] = bse.bandk_dn[k].values[c] - bse.bandk_dn[k].values[v] # eV
			eupHu[i] = ê ⋅ bse.Γdata_dn.upHu[:, i]
		end
		eupHu2 = map(i->abs2(eupHu[i]) / ΔE[i]^2, Base.OneTo(Nvck_uu)) # Å^2
		result_dd = broaden_peaks(broadening, ωgrid, ΔE, eupHu2; scale, normalize = true)
		result_dd .*= coff

		for i in 1:Nexc
			AeupHu[i] += sum(ivck->eigen_vck.vectors[ivck+Nvck_uu, i] * eupHu[ivck], Base.OneTo(Nvck_dd))
		end
	end

	AeupHu2 = map(i -> abs2(AeupHu[i]) / eigen_vck.values[i]^2, Base.OneTo(Nexc))
	result_wieh = broaden_peaks(broadening, ωgrid, eigen_vck.values, AeupHu2; scale, normalize = true)
	result_wieh .*= coff

	if !isnothing(outfolder)
		mkpath(outfolder)
		_write_absorption(joinpath(outfolder, "absorption_noeh_up.dat"), ωgrid, result_uu)
		_write_absorption(joinpath(outfolder, "absorption_noeh_dn.dat"), ωgrid, result_dd)
		_write_absorption(joinpath(outfolder, "absorption_eh.dat"), ωgrid, result_wieh)
	end
	return result_wieh, result_uu, result_dd
end
