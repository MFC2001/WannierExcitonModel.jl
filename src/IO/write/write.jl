
function Base.write(path::AbstractString, data, ::Type{T}; mode = "w", kwargs...) where {T <: FileFormat}
	mkpath(dirname(path))
	return open(path, mode) do io
		write(io, data, T; kwargs...)
	end
end

include("./BAND.jl")
include("./POSCAR.jl")
include("./QEcell.jl")
include("./TB.jl")
include("./wannier/wannier.jl")
