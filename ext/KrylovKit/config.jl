
export eigsolveconfigure!

const eigsolveconfig = Dict{Symbol, Any}(
	:howmany => 1,
	:which => EigSorter(abs; rev = false),
	:verbosity => 0,
	:tol => KrylovDefaults.tol[],
	:krylovdim => KrylovDefaults.krylovdim[],
	:maxiter => 300,
	:orth => KrylovDefaults.orth,
	:issymmetric => false,
	:ishermitian => false,
	:eager => false,
)
"""
	eigsolveconfigure!(; kwards...)

Set parameters for KrylovKit.eigsolve.
Available parameters are:
	`howmany` => number of eigenvalues to compute. Default is 1.
	`which` => the order of eigenvalues. Default is EigSorter(abs; rev = false).
	`verbosity` => verbosity level. Default is 0.
	`tol` => tolerance for convergence. Default is KrylovDefaults.tol.
	`krylovdim` => dimension of the Krylov subspace. Default is KrylovDefaults.krylovdim.
	`maxiter` => maximum number of iterations. Default is 300.
	`orth` => orthogonalization method. Default is KrylovDefaults.orth.
	`issymmetric` => whether the matrix is symmetric. Default is false.
	`ishermitian` => whether the matrix is Hermitian. Default is false.
	`eager` => whether to use eager evaluation. Default is false.
Detailed description of parameters can be found in the [KrylovKit documentation](https://github.com/Jutho/KrylovKit.jl).
"""
function eigsolveconfigure!(; kwards...)
	for (k, v) in kwards
		eigsolveconfig[k] = v
		if k == :howmany
			eigsolveconfig[:krylovdim] = max(2 * v + 10, 30)
		end
	end
end
