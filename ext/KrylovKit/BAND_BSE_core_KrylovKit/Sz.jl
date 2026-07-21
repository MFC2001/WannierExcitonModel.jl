function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSESz, ::Val{true}, ::Val{:Bloch})
	Nq = length(qpoints)
	BSEband_0 = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_p1 = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_n1 = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	E0 = Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}(Vector{Float64}(undef, 0), Matrix{ComplexF64}(undef, 0, 0))

	N0, Np1, Nn1 = bse(Val(:dimension))
	H0 = Matrix{ComplexF64}(undef, N0, N0)
	Hp1 = Matrix{ComplexF64}(undef, Np1, Np1)
	Hn1 = Matrix{ComplexF64}(undef, Nn1, Nn1)

	for (qi, q) in enumerate(qpoints)
		(H_0, H_p1, H_n1) = bse(H0, Hp1, Hn1, q)
		BSEband_0[qi] = N0 > 0 ? _eigsolve_Hmat(H_0) : E0
		BSEband_p1[qi] = Np1 > 0 ? _eigsolve_Hmat(H_p1) : E0
		BSEband_n1[qi] = Nn1 > 0 ? _eigsolve_Hmat(H_n1) : E0
		println("$qi, band of $q calculation end")
		flush(stdout)
	end
	return BSEband_0, BSEband_p1, BSEband_n1
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSESz, ::Val{false}, ::Val{:Bloch})

	(BSEband_0, BSEband_p1, BSEband_n1) = BAND_BSE(KrylovKit_BSEeigenStrategy, qpoints, bse, Val(true), Val(:Bloch))

	BSEband_0 = _eigen2val(BSEband_0)
	BSEband_p1 = _eigen2val(BSEband_p1)
	BSEband_n1 = _eigen2val(BSEband_n1)

	return BSEband_0, BSEband_p1, BSEband_n1
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSESz, ::Val{true}, ::Val{:Periodic})

	e_kR = [cispi(-2 * (k ⋅ R)) for R in bse.unitcell, k in bse.kgrid]

	Nq = length(qpoints)
	BSEband_0 = Vector{BSE_uEigen}(undef, Nq)
	BSEband_p1 = Vector{BSE_uEigen}(undef, Nq)
	BSEband_n1 = Vector{BSE_uEigen}(undef, Nq)

	N0, Np1, Nn1 = bse(Val(:dimension))
	H0 = Matrix{ComplexF64}(undef, N0, N0)
	Hp1 = Matrix{ComplexF64}(undef, Np1, Np1)
	Hn1 = Matrix{ComplexF64}(undef, Nn1, Nn1)

	for (qi, q) in enumerate(qpoints)
		(H_0, H_p1, H_n1) = bse(H0, Hp1, Hn1, q)
		N0 > 0 && (BSEband_0[qi] = BSE_uEigen(bse, q, _eigsolve_Hmat(H_0), e_kR, Val(0)))
		Np1 > 0 && (BSEband_p1[qi] = BSE_uEigen(bse, q, _eigsolve_Hmat(H_p1), e_kR, Val(1)))
		Nn1 > 0 && (BSEband_n1[qi] = BSE_uEigen(bse, q, _eigsolve_Hmat(H_n1), e_kR, Val(-1)))
		println("$qi, band of $q calculation end")
		flush(stdout)
	end

	return BSEband_0, BSEband_p1, BSEband_n1
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSESz, ::Val{false}, ::Val{:Periodic})
	return BAND_BSE(KrylovKit_BSEeigenStrategy, qpoints, bse, Val(false), Val(:Bloch))
end
function BAND_BSE(::Type{KrylovKit_BSEeigenStrategy}, qpoints, bse::BSESz, ::Val{true}, ::Val{:BlochPeriodic})

	e_kR = [cispi(-2 * (k ⋅ R)) for R in bse.unitcell, k in bse.kgrid]

	Nq = length(qpoints)
	BSEband_0 = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_p1 = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_n1 = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	E0 = Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}(Vector{Float64}(undef, 0), Matrix{ComplexF64}(undef, 0, 0))
	BSEband_u_0 = Vector{BSE_uEigen}(undef, Nq)
	BSEband_u_p1 = Vector{BSE_uEigen}(undef, Nq)
	BSEband_u_n1 = Vector{BSE_uEigen}(undef, Nq)

	N0, Np1, Nn1 = bse(Val(:dimension))
	H0 = Matrix{ComplexF64}(undef, N0, N0)
	Hp1 = Matrix{ComplexF64}(undef, Np1, Np1)
	Hn1 = Matrix{ComplexF64}(undef, Nn1, Nn1)

	for (qi, q) in enumerate(qpoints)
		(H_0, H_p1, H_n1) = bse(H0, Hp1, Hn1, q)
		BSEband_0[qi] = N0 > 0 ? _eigsolve_Hmat(H_0) : E0
		BSEband_p1[qi] = Np1 > 0 ? _eigsolve_Hmat(H_p1) : E0
		BSEband_n1[qi] = Nn1 > 0 ? _eigsolve_Hmat(H_n1) : E0
		N0 > 0 && (BSEband_u_0[qi] = BSE_uEigen(bse, q, BSEband_0[qi], e_kR, Val(0)))
		Np1 > 0 && (BSEband_u_p1[qi] = BSE_uEigen(bse, q, BSEband_p1[qi], e_kR, Val(1)))
		Nn1 > 0 && (BSEband_u_n1[qi] = BSE_uEigen(bse, q, BSEband_n1[qi], e_kR, Val(-1)))
		println("$qi, band of $q calculation end")
		flush(stdout)
	end

	return BSEband_0, BSEband_p1, BSEband_n1, BSEband_u_0, BSEband_u_p1, BSEband_u_n1
end
