function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSESU2, ::Val{true}, ::Val{:Bloch})
	Nq = length(qpoints)
	BSEband_t = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_s = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)

	N = length(bse.vckmap)
	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)
		(H_t, H_s) = bse(Htriplet, Hsinglet, q)
		BSEband_t[qi] = eigen!(H_t)
		BSEband_s[qi] = eigen!(H_s)
	end
	return BSEband_t, BSEband_s
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSESU2, ::Val{false}, ::Val{:Bloch})
	Nq = length(qpoints)
	N = length(bse.vckmap)

	BSEband_t = Matrix{Float64}(undef, N, Nq)
	BSEband_s = Matrix{Float64}(undef, N, Nq)

	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)
		(H_t, H_s) = bse(Htriplet, Hsinglet, q)
		BSEband_t[:, qi] = eigvals!(H_t)
		BSEband_s[:, qi] = eigvals!(H_s)
	end
	return BSEband_t, BSEband_s
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSESU2, ::Val{true}, ::Val{:Periodic})

	e_kR = [cis(-2π * (k ⋅ R)) for R in bse.unitcell, k in bse.kgrid]

	Nq = length(qpoints)
	BSEband_t = Vector{BSE_uEigen}(undef, Nq)
	BSEband_s = Vector{BSE_uEigen}(undef, Nq)

	N = length(bse.vckmap)
	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)
		(H_t, H_s) = bse(Htriplet, Hsinglet, q)
		BSEband_t[qi] = BSE_uEigen(bse, q, eigen!(H_t), e_kR)
		BSEband_s[qi] = BSE_uEigen(bse, q, eigen!(H_s), e_kR)
	end

	return BSEband_t, BSEband_s
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSESU2, ::Val{false}, ::Val{:Periodic})
	return BAND_BSE(LinearAlgebra_BSEeigenStrategy, qpoints, bse, Val(false), Val(:Bloch))
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSESU2, ::Val{true}, ::Val{:BlochPeriodic})

	e_kR = [cis(-2π * (k ⋅ R)) for R in bse.unitcell, k in bse.kgrid]

	Nq = length(qpoints)
	BSEband_t = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_s = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_u_t = Vector{BSE_uEigen}(undef, Nq)
	BSEband_u_s = Vector{BSE_uEigen}(undef, Nq)

	N = length(bse.vckmap)
	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)
		(H_t, H_s) = bse(Htriplet, Hsinglet, q)
		BSEband_t[qi] = eigen!(H_t)
		BSEband_s[qi] = eigen!(H_s)
		BSEband_u_t[qi] = BSE_uEigen(bse, q, BSEband_t[qi], e_kR)
		BSEband_u_s[qi] = BSE_uEigen(bse, q, BSEband_s[qi], e_kR)
	end

	return BSEband_t, BSEband_u_t, BSEband_s, BSEband_u_s
end
