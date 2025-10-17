function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEcluster_SU2, ::Val{true})
	BSEband_t = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, 1)
	BSEband_s = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, 1)

	(H_t, H_s) = bse()
	BSEband_t[1] = eigen!(H_t)
	BSEband_s[1] = eigen!(H_s)

	return BSEband_t, BSEband_s
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEcluster_SU2, ::Val{false})
	N = length(bse.vcmap)

	BSEband_t = Matrix{Float64}(undef, N, 1)
	BSEband_s = Matrix{Float64}(undef, N, 1)

	(H_t, H_s) = bse()
	BSEband_t[:, 1] = eigvals!(H_t)
	BSEband_s[:, 1] = eigvals!(H_s)

	return BSEband_t, BSEband_s
end
