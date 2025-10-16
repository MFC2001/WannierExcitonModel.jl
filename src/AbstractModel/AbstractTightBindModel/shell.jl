
function Kline(nk::Integer, TB::AbstractTightBindModel)
	return Kline(nk, TB.lattice, TB.period)
end
function Kline(point::AbstractVector{<:Pair}, nk::Integer, TB::AbstractTightBindModel)
	return Kline(point, nk, TB.lattice)
end

Cell(TB::TightBindModel) = Cell{ReducedCoordinates{Float64}}(TB.lattice, TB.atom_location, TB.atom_name, eachindex(TB.atom_location), TB.period)
HR(TB::TightBindModel) = HR(TB.H)
function ORBITAL(TB::AbstractTightBindModel)
	return ORBITAL(
		map(x -> TB.lattice * x, TB.orb_location);
		name = TB.orb_name,
		atom_location = map(x -> TB.lattice * x, TB.atom_location),
		atom_name = TB.atom_name,
	)
end

"""
	UwithLR(U, TB::TightBindModel, rcut::Real; kwargs...) -> UwithLR

A shortcut to create a `UwithLR`.
This long-range correction is achieved by using Gauss potential.

- `U`: the direct Coulomb potential between the wannier bases within a limited distance, can be a `path::String`, `hr::HR` or `rh::AbstractReciprocalHoppings`.
kwargs:
- `ϵ::Real`: dielectric constant;
- `αrcut::Real`: α × rcut, default is 4.5;
- `δ:Real`: the cutoff value of reciprocal potential.
"""
function UwithLR(U::AbstractReciprocalHoppings, TB::AbstractTightBindModel, rcut::Real; kwards...)
	return UwithLR(U, TB.lattice, TB.orb_location, rcut; period = TB.period, kwards...)
end
"""
	MirrorCorrection(U::HR, TB::TightBindModel, rcut::Real; kwargs...) -> HR

Try to apply a mirror correction to `U`.
This mirror correction is achieved by using Gauss potential.

- `rcut::Real`: should be less than a real radius value, out of which the U approximate to 1/r, default is judged by `U` and `TB`;
kwargs:
- `kgrid::MonkhorstPack`: should be the kgrid used to calculate `U`, default is judged by `U`;
- `ϵ::Real`: dielectric constant;
- `αrcut::Real`: α × rcut, default is 4.5;
- `δ:Real`: the cutoff value of reciprocal potential.
"""
function MirrorCorrection(U::HR, TB::AbstractTightBindModel, rcut::Union{Nothing, Real} = nothing; kwargs...)
	return MirrorCorrection(U, TB.lattice, TB.orb_location, rcut; period = TB.period, kwargs...)
end
