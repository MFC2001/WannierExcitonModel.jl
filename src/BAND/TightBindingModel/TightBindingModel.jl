
include("./core.jl")

"""
	BAND(kpoints, TB::AbstractTightBindModel; vector::Bool = false, wfctype::Symbol = :Bloch)

This method can calculate the band structure of a tight-binding model.
- `kpoints` can be a `AbstractBrillouinZone` object, a vector of kpoints or a single kpoint, where kpoint should be a vector of 3 real numbers.
- `vector` can be set as true, which decides whether to output wavefunction.
- `wfctype` can be set as :Periodic, which decides the type of wavefunction.
"""
function BAND(kpoints::AbstractVector{<:ReducedCoordinates}, TB::AbstractTightBindModel; vector::Bool = false, wfctype::Symbol = :Bloch)
	return BAND_TB(kpoints, TB, Val(wfctype), vector)
end
function BAND_TB(kpoints::AbstractVector{<:ReducedCoordinates}, TB::AbstractTightBindModel, ::Val{:Bloch}, vector::Bool)
	return BAND_hrh(kpoints, TB.H, Val(vector))
end
function BAND_TB(kpoints::AbstractVector{<:ReducedCoordinates}, TB::AbstractTightBindModel, ::Val{:Periodic}, vector::Bool)
	return BAND_hrh(kpoints, TB.H, TB.orb_location, Val(vector))
end
"""
	BAND(kpoints::AbstractVector{<:ReducedCoordinates}, hr::HR; vector::Bool = false)
	BAND(kpoints::AbstractVector{<:ReducedCoordinates}, hr::HR, orblocat::AbstractVector{<:ReducedCoordinates}; vector::Bool = false)

The second function will return the periodic part of the wavefunction if `vector` = true.
"""
function BAND(kpoints::AbstractVector{<:ReducedCoordinates}, hr::HR; vector::Bool = false)
	return BAND_hrh(kpoints, HermitianReciprocalHoppings(hr), Val(vector))
end
function BAND(kpoints::AbstractVector{<:ReducedCoordinates}, hr::HR, orblocat::AbstractVector{<:ReducedCoordinates}; vector::Bool = false)
	return BAND_hrh(kpoints, HermitianReciprocalHoppings(hr), orblocat, Val(vector))
end
"""
	BAND(kpoints::AbstractVector{<:ReducedCoordinates}, rh::HermitianReciprocalHoppings; vector::Bool = false)
	BAND(kpoints::AbstractVector{<:ReducedCoordinates}, rh::HermitianReciprocalHoppings, orblocat::AbstractVector{<:ReducedCoordinates}; vector::Bool = false)

This two function will return the bloch wavefunction if `vector` = true.
"""
function BAND(kpoints::AbstractVector{<:ReducedCoordinates}, rh::HermitianReciprocalHoppings; vector::Bool = false)
	return BAND_hrh(kpoints, rh, Val(vector))
end
function BAND(kpoints::AbstractVector{<:ReducedCoordinates}, rh::HermitianReciprocalHoppings, orblocat::AbstractVector{<:ReducedCoordinates}; vector::Bool = false)
	return BAND_hrh(kpoints, rh, orblocat, Val(vector))
end

"""
	ExtractU(kpoints, band, TB::AbstractTightBindModel) -> band_u
	ExtractU(kpoints, band, orblocat::AbstractVector{<:ReducedCoordinates}) -> band_u

This function can extract the periodic part of the wavefunction from the bloch wavefunction.
- `kpoints` can be a `AbstractBrillouinZone` object, a vector of kpoints or a single kpoint, where kpoint should be a vector of 3 real numbers.
- `band` can be a vector of `Eigen` object or a single `Eigen` object.
- `band_u` is a vector of `Eigen` object or a single `Eigen` object depending on the input `band`.
"""
function ExtractU(kpoints, band, TB::AbstractTightBindModel)
	return ExtractU(kpoints, band, TB.orb_location)
end
function ExtractU(kpoints::AbstractVector{<:ReducedCoordinates}, band::AbstractVector{<:Eigen}, orblocat::AbstractVector{<:ReducedCoordinates})
	length(kpoints) == length(band) || error("Mismatched kpoints and band.")
	band_u = similar(band)
	for k in eachindex(kpoints)
		band_u[k] = ExtractU(kpoints[k], band[k], orblocat)
	end
	return band_u
end
function ExtractU(kpoint::ReducedCoordinates, band::Eigen, orblocat::AbstractVector{<:ReducedCoordinates})
	return Eigen(copy(band.values), ExtractU_TB(kpoint, band.vectors, orblocat))
end
function ExtractU_TB(kpoint, band, orblocat)
	ekτ = map(τ -> cis(-2π * (kpoint ⋅ τ)), orblocat)
	return ekτ .* band
end
