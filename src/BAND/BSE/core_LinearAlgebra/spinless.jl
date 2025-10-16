function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEspinless, ::Val{true}, ::Val{:Bloch}; ηt, ηs)
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
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEspinless, ::Val{false}, ::Val{:Bloch}; ηt, ηs)
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
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEspinless, ::Val{true}, ::Val{:Periodic}; ηt, ηs)

	ijRvck = _uijR_ψvck(bse, ηt, ηs)

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
		BMt, BMs = ijRvck(bse.bandk, bse.bandkq, q)
		BSEband_t[qi] = Eigen(BSEband_t[qi].values, BMt * BSEband_t[qi].vectors)
		BSEband_s[qi] = Eigen(BSEband_s[qi].values, BMs * BSEband_s[qi].vectors)
	end

	return BSEband_t, BSEband_s
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEspinless, ::Val{false}, ::Val{:Periodic}; ηt, ηs)
	return BAND_BSE(LinearAlgebra_BSEeigenStrategy, qpoints, bse, Val(false), Val(:Bloch); ηt, ηs)
end
function BAND_BSE(::Type{LinearAlgebra_BSEeigenStrategy}, qpoints, bse::BSEspinless, ::Val{true}, ::Val{:BlochPeriodic}; ηt, ηs)

	ijRvck = _uijR_ψvck(bse, ηt, ηs)

	Nq = length(qpoints)
	BSEband_t = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_s = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_u_t = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)
	BSEband_u_s = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, Nq)

	N = length(bse.vckmap)
	Htriplet = Matrix{ComplexF64}(undef, N, N)
	Hsinglet = Matrix{ComplexF64}(undef, N, N)

	for (qi, q) in enumerate(qpoints)
		(H_t, H_s) = bse(Htriplet, Hsinglet, q)
		BSEband_t[qi] = eigen!(H_t)
		BSEband_s[qi] = eigen!(H_s)
		BMt, BMs = ijRvck(bse.bandk, bse.bandkq, q)
		BSEband_u_t[qi] = Eigen(copy(BSEband_t[qi].values), BMt * BSEband_t[qi].vectors)
		BSEband_u_s[qi] = Eigen(copy(BSEband_s[qi].values), BMs * BSEband_s[qi].vectors)
	end

	return BSEband_t, BSEband_u_t, BSEband_s, BSEband_u_s
end
