function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEcluster_general, ::Val{true})
	BSEband = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, 1)

	H_ = bse()
	BSEband[1] = _eigsolve_Hmat(H_)

	return BSEband
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEcluster_general, ::Val{false})
	BSEband = BAND_BSE(KrylovKit_BSEeigenStrategy, qpoints, bse, Val(true))
	BSEband = _eigen2vals(BSEband)
	return BSEband
end
