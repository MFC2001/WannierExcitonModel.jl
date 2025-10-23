function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{true}, ::Val{:Bloch})
	Nq = length(qpoints)
	BSEband = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)

	N = length(bse.vckmap)
	H = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)
		H_ = bse(H, q)
		BSEband[qi] = _eigsolve_Hmat(H_)
	end
	return BSEband
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{false}, ::Val{:Bloch})
	BSEband = BAND_BSE(KrylovKit_BSEeigenStrategy, qpoints, bse, Val(true), Val(:Bloch))
	BSEband = _eigen2vals(BSEband)
	return BSEband
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{true}, ::Val{:Periodic})

	e_kR = [cis(-2π * (k ⋅ R)) for R in bse.unitcell, k in bse.kgrid]

	Nq = length(qpoints)
	BSEband = Vector{BSE_uEigen}(undef, Nq)

	N = length(bse.vckmap)
	H = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)
		H_ = bse(H, q)
		BSEband[qi] = BSE_uEigen(bse, q, _eigsolve_Hmat(H_), e_kR)
	end

	return BSEband
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{false}, ::Val{:Periodic})
	return BAND_BSE(KrylovKit_BSEeigenStrategy, qpoints, bse, Val(false), Val(:Bloch))
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSEgeneral, ::Val{true}, ::Val{:BlochPeriodic})

	e_kR = [cis(-2π * (k ⋅ R)) for R in bse.unitcell, k in bse.kgrid]

	Nq = length(qpoints)
	BSEband = Vector{BSE_uEigen}(undef, Nq)
	BSEband_u = Vector{BSE_uEigen}(undef, Nq)

	N = length(bse.vckmap)
	H = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)
		H_ = bse(H, q)
		BSEband[qi] = _eigsolve_Hmat(H_)
		BSEband_u[qi] = BSE_uEigen(bse, q, BSEband[qi], e_kR)
	end

	return BSEband, BSEband_u
end
