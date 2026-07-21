
function WannierInterActionBase.Kline(nk::Integer, TB::AbstractTightBindModel)
	return Kline(nk, TB.lattice, TB.period)
end
function WannierInterActionBase.Kline(point::AbstractVector{<:Pair}, nk::Integer, TB::AbstractTightBindModel)
	return Kline(point, nk, TB.lattice)
end

WannierInterActionBase.Cell(TB::TightBindModel) = Cell{ReducedCoordinates{Float64}}(TB.lattice, TB.atom_location, TB.atom_name, eachindex(TB.atom_location), TB.period)
WannierInterActionBase.HR(TB::TightBindModel) = HR(TB.H)
function WannierInterActionBase.wannier90_centres(TB::AbstractTightBindModel)
	return wannier90_centres(
		map(x -> TB.lattice * x, TB.orb_location);
		name = TB.orb_name,
		atom_location = map(x -> TB.lattice * x, TB.atom_location),
		atom_name = TB.atom_name,
	)
end
function WannierInterActionBase.translate(TB::TightBindModel, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)
	new_orblocat = copy(TB.orb_location)
	for (iw, R) in centres
		new_orblocat[iw] = new_orblocat[iw] - ReducedCoordinates{Int}(R)
	end
	new_H = translate(TB.H, centres...)
	return TightBindModel(TB.lattice, copy(TB.atom_name), copy(TB.atom_location), copy(TB.orb_name), new_orblocat,
		new_H, TB.period)
end
function WannierInterActionBase.translate!(TB::TightBindModel, centres::Pair{<:Integer, <:AbstractVector{<:Integer}}...)
	for (iw, R) in centres
		TB.orb_location[iw] = TB.orb_location[iw] - ReducedCoordinates{Int}(R)
	end
	translate!(TB.H, centres...)
	return TB
end

# """
# 	UwithLR(U, TB::TightBindModel, rcut::Real; kwargs...) -> UwithLR

# A shortcut to create a `UwithLR`.
# This long-range correction is achieved by using Gauss potential.

# - `U`: the direct Coulomb potential between the wannier bases within a limited distance, can be a `path::String`, `hr::HR` or `rh::AbstractReciprocalHoppings`.
# kwargs:
# - `ϵ::Real`: dielectric constant;
# - `αrcut::Real`: α × rcut, default is 4.5;
# - `δ:Real`: the cutoff value of reciprocal potential.
# """
# function WannierInterActionBase.UwithLR(U::AbstractReciprocalHoppings, TB::AbstractTightBindModel, rcut::Real; kwards...)
# 	return WannierInterActionBase.UwithLR(U, TB.lattice, TB.orb_location, rcut; period = TB.period, kwards...)
# end
# """
# 	MirrorCorrection(U::HR, TB::TightBindModel, rcut::Real; kwargs...) -> HR

# Try to apply a mirror correction to `U`.
# This mirror correction is achieved by using Gauss potential.

# - `rcut::Real`: should be less than a real radius value, out of which the U approximate to 1/r, default is judged by `U` and `TB`;
# kwargs:
# - `kgrid::MonkhorstPack`: should be the kgrid used to calculate `U`, default is judged by `U`;
# - `ϵ::Real`: dielectric constant;
# - `αrcut::Real`: α × rcut, default is 4.5;
# - `δ:Real`: the cutoff value of reciprocal potential.
# """
# function WannierInterActionBase.MirrorCorrection(U::HR, TB::AbstractTightBindModel, rcut::Union{Nothing, Real} = nothing; kwargs...)
# 	return MirrorCorrection(U, TB.lattice, TB.orb_location, rcut; period = TB.period, kwargs...)
# end
