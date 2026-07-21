export absorption

include("./SU2_general.jl")
include("./Sz.jl")

function _write_absorption(file, ωgrid, data)
	io = open(file, "w")
	println(io, " # Column 1: omega")
	println(io, " # Column 2: eps2(omega)")
	for (ω, ϵ₂) in zip(ωgrid, data)
		@printf(io, "%16.9f %16.9f\n", ω, ϵ₂)
	end
	close(io)
	return nothing
end
