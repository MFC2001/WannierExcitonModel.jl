function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEcluster_SU2, ::Val{true})

	BSEband_t = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, 1)
	BSEband_s = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, 1)

	(H_t, H_s) = bse()
	BSEband_t[1] = _eigsolve_Hmat(H_t)
	BSEband_s[1] = _eigsolve_Hmat(H_s)

	return BSEband_t, BSEband_s
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEcluster_SU2, ::Val{false})

	(BSEband_t, BSEband_s) = BAND_BSE(KrylovKit_BSEeigenStrategy, qpoints, bse, Val(true))

	BSEband_t = _eigen2vals(BSEband_t)
	BSEband_s = _eigen2vals(BSEband_s)

	return BSEband_t, BSEband_s
end
