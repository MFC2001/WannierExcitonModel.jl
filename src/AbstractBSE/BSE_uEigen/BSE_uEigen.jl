export BSE_uEigen
"""
	abstract type BSE_uEigen end

Used to calculate the periodic part of excitonic Bloch function.

	(u::BSE_uEigen)(λ::Real, bandindex::AbstractVector{<:Integer}=1:nband) -> Eigen

Input the weight `λ` of exciton center coordinates, return the numerical vector of periodic part.

```julia
julia> u_λ = u(λ); u_λ isa Eigen

true
```
"""
abstract type BSE_uEigen end

include("./wiTDA/wiTDA.jl")
include("./woTDA/woTDA.jl")
