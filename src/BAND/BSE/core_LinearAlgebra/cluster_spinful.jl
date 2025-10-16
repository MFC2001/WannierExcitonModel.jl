function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEcluster_spinful, ::Val{true})
	BSEband = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, 1)

	H_ = bse()
	BSEband[1] = eigen!(H_)

	return BSEband
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEcluster_spinful, ::Val{false})

	N = length(bse.vckmap)

	BSEband = Matrix{Float64}(undef, N, 1)

	H_ = bse()
	BSEband_t[:, 1] = eigvals!(H_)

	return BSEband
end
