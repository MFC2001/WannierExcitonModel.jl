export FileFormat
export BAND, POSCAR

"""
	read(filename::AbstractString, ::Type{FileFormat}; kwargs...)
	read(io::IO, ::Type{BAND}; kwargs...)
"""
abstract type FileFormat end
struct BAND <: FileFormat end
struct POSCAR <: FileFormat end

include("./HR.jl")
include("./ORBITAL.jl")
