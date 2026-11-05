export BSE, numatom, allocate!, release!, isallocated

"""
	_BSE_buffer

Abstract type of buffer array used to construct BSE Hamiltonian.
"""
abstract type _BSE_buffer end
include("./wiTDA/wiTDA.jl")
include("./woTDA/woTDA.jl")

numatom(bse::AbstractBSE) = length(bse.atom_location)
(bse::AbstractBSE)(sym::Symbol, args...; kwargs...) = bse(Val(sym), args...; kwargs...)
(bse::AbstractBSE)(q::AbstractVector{<:Real}; kwargs...) = bse(ReducedCoordinates(q); kwargs...)
isallocated(bse::AbstractBSE) = bse.buffer.isallocated

"""
Spin degeneracy situation:

	BSE(kgrid, TB, Kᵈ, Kˣ; v::AbstractVector{<:Integer}, c::AbstractVector{<:Integer}, scissor::Real=0)

- `v` and `c` are the indices of valence and conduction bands, respectively.
- `scissor` is the scissor correction to the conduction bands, default is 0.

Noncollinear situation:

	BSE(kgrid, TB, Kᵈ, Kˣ; channel::Integer, v::AbstractVector{<:Integer}, c::AbstractVector{<:Integer}, scissor::Real=0)

- `channel` have to be set as 1.

Spin-polarized situation:

    BSE(kgrid, TB_up, TB_dn, Kᵈ, Kˣ; v_up::AbstractVector{<:Integer}=Int[], c_up::AbstractVector{<:Integer}=Int[],
    v_dn::AbstractVector{<:Integer}=Int[], c_dn::AbstractVector{<:Integer}=Int[],
    scissor::Real=0, scissor_up::Real=scissor, scissor_dn::Real=scissor,
    scissor_c_up::Real=scissor_up, scissor_c_dn::Real=scissor_dn, scissor_v_up::Real=0, scissor_v_dn::Real=0) -> AbstractBSE

- `v_up`, `c_up`, `v_dn`, `c_dn` are the indices of valence and conduction bands for spin-up and spin-down channels, respectively.
- `scissor_c_up` is the energy shift of the conduction bands for spin-up channel, the others are similar as well.

```julia
julia> bse isa AbstractBSE
true

julia> bse([0.1, 0, 0]) isa Hermitian
true
```
"""
function BSE(kgrid::Union{MonkhorstPack,RedKgrid,AbstractVector{<:Integer},NTuple{3,<:Integer}},
    args...; TDA::Bool=true, kwargs...)
    if kgrid isa AbstractVector || kgrid isa NTuple
        kgrid = RedKgrid(MonkhorstPack(kgrid))
    elseif kgrid isa MonkhorstPack
        kgrid = RedKgrid(kgrid)
    end
    @assert length(kgrid) > 0 "Empty `kgrid` is unaccepted!"
    if TDA
        return BSE_wiTDA(kgrid, args...; kwargs...)
    else
        return BSE_woTDA(kgrid, args...; kwargs...)
    end
end

#From Rohlfing and Louie, we need to control the arbitrary phase of electron wavefunction.
# Acctually, the random phase of electron wavefunction does not change the excitonic result;
# But a definite phase is helpful to get the exact electron wavefunction anytime.

function _BSE_preprocess_Γ(q, isΓ, period)
    # Special treatment is required only when q is equivalent to Γ.
    if iszero(q) || all(isinteger, q)
        np = count(period)
        if np ≠ 3
            # only in 3D case, we need to add a term which depends on the direction of q.
            return ReducedCoordinates(0, 0, 0), true
        else
            # 3D case
            if isΓ
                if iszero(q)
                    @info "A zero vector `q` is unacceptable for Γ calculation in 3D system."
                else
                    # q is not zero, then q is direction
                    return ReducedCoordinates(q), true
                end
            else
                @info "A zero vector or an integer vector `q` is unacceptable for non-Γ calculation in 3D system, set `isΓ` as true."
            end
            @info "set ̂q as [0, 0, 1]."
            return ReducedCoordinates(0, 0, 1), true
        end
    else
        return q, isΓ
    end
end
