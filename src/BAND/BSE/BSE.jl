
abstract type BSEeigenStrategy end
struct LinearAlgebra_BSEeigenStrategy <: BSEeigenStrategy end
struct KrylovKit_BSEeigenStrategy <: BSEeigenStrategy end

const _BSE_eigen_STRATEGY = Ref{Type{<:BSEeigenStrategy}}(LinearAlgebra_BSEeigenStrategy)

"""
Set the eigenvalue solving strategy for BSE calculations.
"""
set_BSE_eigen_strategy!(::Type{S}) where {S <: BSEeigenStrategy} = _BSE_eigen_STRATEGY[] = S

include("./core_LinearAlgebra/core_LinearAlgebra.jl")
include("./ExtractU/ExtractU.jl")


"""
	BAND(qpoints, bse::AbstractBSE; vector::Bool = false, wfctype::Symbol = :Bloch, 
		η::Real = 1 // 2, ηt::Real = η, ηs::Real = η)

This function can calculate the band structure of a BSE model.
- `qpoints` can be a `AbstractBrillouinZone` object, a vector of qpoints or a single qpoint, where qpoint should be a vector of 3 real numbers.
- `vector` can be set as true, which decides whether to output wavefunction.
- `wfctype` can be set as :Periodic, which decides the type of wavefunction.
- `η` is the parameter of excitonic periodic wavefunction, also determines the excitonic center position. For spinless BSE model, `ηt` and `ηs` can be set separately.

We provide two strategies to solve the BSE Hamiltonian eigenvalue problem: 

	LinearAlgebra_BSEeigenStrategy <: BSEeigenStrategy
	KrylovKit_BSEeigenStrategy <: BSEeigenStrategy

- `LinearAlgebra_BSEeigenStrategy`: use LinearAlgebra to solve eigenvalue equation;
- `KrylovKit_BSEeigenStrategy`: use KrylovKit to solve eigenvalue equation.
The default strategy is `LinearAlgebra_BSEeigenStrategy`, which uses the built-in `eigen` function in Julia.
When `KrylovKit` package is used, the default strategy is set as `KrylovKit_BSEeigenStrategy` automatically.

```julia
julia> using KrylovKit; using ExcitonicWannierModel;
```

You can switch the strategy by calling `set_BSE_eigen_strategy!(::Type{BSEeigenStrategy})` manually.

When `KrylovKit_BSEeigenStrategy` is used, you can set parameters of `KrylovKit.eigsolve` by calling `eigsolveconfigure!(; kwards...)`.
For excitonic Hamiltonian here, the relevant kwargs are: `howmany`, `which`, `verbosity`, `tol`, `krylovdim`, `maxiter`, `orth`.
"""
function BAND(qpoints::AbstractVector{<:ReducedCoordinates}, bse::BSEspinless;
	vector::Bool = false, wfctype::Symbol = :Bloch, η::Real = 1 // 2, ηt::Real = η, ηs::Real = η)
	return BAND_BSE(_BSE_eigen_STRATEGY[], qpoints, bse, Val(vector), Val(wfctype); ηt, ηs)
end
function BAND(qpoints::AbstractVector{<:ReducedCoordinates}, bse::BSEspinful;
	vector::Bool = false, wfctype::Symbol = :Bloch, η::Real = 1 // 2)
	return BAND_BSE(_BSE_eigen_STRATEGY[], qpoints, bse, Val(vector), Val(wfctype); η)
end
function BAND(qpoints::AbstractVector{<:ReducedCoordinates}, bse::BSEcluster_spinless; vector::Bool = false)
	return BAND_BSE(_BSE_eigen_STRATEGY[], qpoints, bse, Val(vector))
end
function BAND(qpoints::AbstractVector{<:ReducedCoordinates}, bse::BSEcluster_spinful; vector::Bool = false)
	return BAND_BSE(_BSE_eigen_STRATEGY[], qpoints, bse, Val(vector))
end
