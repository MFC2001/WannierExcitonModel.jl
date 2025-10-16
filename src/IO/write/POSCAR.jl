
function Base.write(io::IO, cell::Cell, ::Type{POSCAR}; comment = "", format = "Cartesian")
	now = Dates.format(Dates.now(), "yyyy/m/d H:M:S")
	comment = comment * " " * now

	println(io, comment)
	println(io, "1.0")

	for i in 1:3
		@printf(io, "%26.16f %26.16f %26.16f\n", cell.lattice[:, i]...)
	end

	elem_name = atomnames(cell)
	for name in elem_name
		@printf(io, " %5s", name)
	end
	println(io, "")

	for name in elem_name
		@printf(io, " %5u", count(cell.name .== name))
	end
	println(io, "")

	if format[1] ∈ ['C', 'c']
		write(io, "Cartesian\n")
		for r in cell.location
			@printf(io, "%23.16f %23.16f %23.16f\n", cell.lattice * r...)
		end
	elseif format[1] ∈ ['D', 'd']
		write(io, "Direct\n")
		for r in cell.location
			@printf(io, "%23.16f %23.16f %23.16f\n", r...)
		end
	end

	return nothing
end

Base.write(io::IO, TB::AbstractTightBindModel, ::Type{POSCAR}; kwargs...) =
	write(io, Cell(TB), POSCAR; kwargs...)
