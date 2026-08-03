export BAND, ExtractU, ExpectationValue

include("./TightBindingModel/TightBindingModel.jl")
include("./BSE/BSE.jl")

function BAND(kgrid::MonkhorstPack, args...; kwargs...)
	return BAND(RedKgrid(kgrid), args...; kwargs...)
end
function BAND(kpoints::AbstractBrillouinZone, args...; kwargs...)
	return BAND(kpoints.kdirect, args...; kwargs...)
end
function BAND(kpoints::AbstractVector{<:AbstractVector{<:Real}}, args...; kwargs...)
	return BAND(ReducedCoordinates.(kpoints), args...; kwargs...)
end
function BAND(kpoint::AbstractVector{<:Real}, args...; kwargs...)
	return BAND([ReducedCoordinates(kpoint)], args...; kwargs...)
end

function ExtractU(kgrid::MonkhorstPack, args...; kwargs...)
	return ExtractU(RedKgrid(kgrid), args...; kwargs...)
end
function ExtractU(kpoints::AbstractBrillouinZone, args...; kwargs...)
	return ExtractU(kpoints.kdirect, args...; kwargs...)
end
function ExtractU(kpoints::AbstractVector{<:AbstractVector{<:Real}}, args...; kwargs...)
	return ExtractU(ReducedCoordinates.(kpoints), args...; kwargs...)
end
function ExtractU(kpoint::AbstractVector{<:Real}, args...; kwargs...)
	return ExtractU(ReducedCoordinates(kpoint), args...; kwargs...)
end

"""
	ExpectationValue(operator::AbstractMatrix{<:Number}, band::AbstractVector{<:Eigen}) -> Vector{Vector{Float64}}
	ExpectationValue(operator::AbstractMatrix{<:Number}, band::Eigen) -> Vector{Float64}

Calculate the expectation value of an operator for each eigenstate.
"""
function ExpectationValue(operator::AbstractMatrix{<:Number}, band::AbstractVector{<:Eigen})
	size(band[1].vectors, 1) == size(operator, 1) == size(operator, 2) || error("Mismatched operator and band.")
	expect_val = similar(band, Vector{Float64})
	for i in eachindex(band)
		expect_val[i] = map(eachcol(band[i].vectors)) do ψ
			real(ψ ⋅ (operator * ψ))
		end
	end
	return expect_val
end
function ExpectationValue(operator::AbstractMatrix{<:Number}, band::Eigen)
	size(band.vectors, 1) == size(operator, 1) == size(operator, 2) || error("Mismatched operator and band.")
	expect_val = map(eachcol(band.vectors)) do ψ
		real(ψ ⋅ (operator * ψ))
	end
	return expect_val
end
