function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{true}, ::Val{:Bloch}; η)
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
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{false}, ::Val{:Bloch}; η)
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
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{true}, ::Val{:Periodic}; η)

	ijRvck = _uijR_ψvck(bse, η)

	Nq = length(qpoints)
	BSEband = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)

	N = length(bse.vckmap)
	H = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)
		H_ = bse(H, q)
		BSEband[qi] = eigen!(H_)
		BM = ijRvck(bse.bandk, bse.bandkq, q)
		BSEband[qi] = Eigen(BSEband[qi].values, BM * BSEband[qi].vectors)
	end

	return BSEband
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{false}, ::Val{:Periodic}; η)
	return BAND_BSE(LinearAlgebra_BSEeigenStrategy, qpoints, bse, Val(false), Val(:Bloch); η)
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{true}, ::Val{:BlochPeriodic}; η)

	ijRvck = _uijR_ψvck(bse, η)

	Nq = length(qpoints)
	BSEband = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_u = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)

	N = length(bse.vckmap)
	H = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)
		H_ = bse(H, q)
		BSEband[qi] = eigen!(H_)
		BM = ijRvck(bse.bandk, bse.bandkq, q)
		BSEband_u[qi] = Eigen(copy(BSEband[qi].values), BM * BSEband[qi].vectors)
	end

	return BSEband, BSEband_u
end
