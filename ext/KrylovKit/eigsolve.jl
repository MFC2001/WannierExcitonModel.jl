
function _eigsolve_Hmat(H)

	(vals, vecs, info) = KrylovKit.eigsolve(H, eigsolveconfig[:howmany], eigsolveconfig[:which];
		verbosity = eigsolveconfig[:verbosity],
		tol = eigsolveconfig[:tol],
		krylovdim = min(eigsolveconfig[:krylovdim], size(H, 1)),
		maxiter = eigsolveconfig[:maxiter],
		orth = eigsolveconfig[:orth],
		ishermitian = true,
	)

	if info.converged < eigsolveconfig[:howmany]
		@warn "KrylovKit.eigsolve don't converge!"
	elseif info.converged > eigsolveconfig[:howmany]
		vals = vals[1:eigsolveconfig[:howmany]]
		vecs = vecs[1:eigsolveconfig[:howmany]]
	end

	I = sortperm(vals)
	vals = vals[I]

	vecs_mat = Matrix{eltype(eltype(vecs))}(undef, size(H, 1), eigsolveconfig[:howmany])
	for i in eachindex(I)
		vecs_mat[:, i] = vecs[I[i]]
	end

	return Eigen(vals, vecs_mat)
end
