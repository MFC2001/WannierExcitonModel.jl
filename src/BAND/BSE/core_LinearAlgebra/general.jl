function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{true}, ::Val{:Bloch})
	Nq = length(qpoints)
	BSEband = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)

	N = length(bse.vckmap)
	H = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)
		H_ = bse(H, q)
		BSEband[qi] = eigen!(H_)
	end
	return BSEband
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{false}, ::Val{:Bloch})
	Nq = length(qpoints)
	N = length(bse.vckmap)

	BSEband = Matrix{Float64}(undef, N, Nq)

	H = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)
		H_ = bse(H, q)
		BSEband_t[:, qi] = eigvals!(H_)
	end
	return BSEband
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{true}, ::Val{:Periodic})

	e_kR = [cis(-2π * (k ⋅ R)) for R in bse.unitcell, k in bse.kgrid]

	Nq = length(qpoints)
	BSEband = Vector{BSE_uEigen}(undef, Nq)

	N = length(bse.vckmap)
	H = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)
		H_ = bse(H, q)
		BSEband[qi] = BSE_uEigen(bse, q, eigen!(H_), e_kR)
	end

	return BSEband
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{false}, ::Val{:Periodic}; η)
	return BAND_BSE(LinearAlgebra_BSEeigenStrategy, qpoints, bse, Val(false), Val(:Bloch); η)
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{true}, ::Val{:BlochPeriodic}; η)

	e_kR = [cis(-2π * (k ⋅ R)) for R in bse.unitcell, k in bse.kgrid]

	Nq = length(qpoints)
	BSEband = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_u = Vector{BSE_uEigen}(undef, Nq)

	N = length(bse.vckmap)
	H = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)
		H_ = bse(H, q)
		BSEband[qi] = eigen!(H_)
		BSEband_u[qi] = BSE_uEigen(bse, q, BSEband[qi], e_kR)
	end

	return BSEband, BSEband_u
end
