
function Base.read(path::AbstractString, ::Type{T}; kwargs...) where {T <: FileFormat}
	return open(path, "r") do io
		read(io, T; kwargs...)
	end
end

include("./BAND.jl")
include("./POSCAR.jl")
include("./RESPACKUJ.jl")
include("./dat_eigenvalue.jl")
include("./wannier/wannier.jl")
