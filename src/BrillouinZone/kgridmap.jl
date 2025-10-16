export kgridmap
"""
	kgridmap(aimkdirects, kdirects, mode = +) -> Matrix{Int}

Try to get a map matrix between aimkdirects and kdirects.

```julia
julia> M = kgridmap(aimkdirects, kdirects, +);

julia> aimkdirects[M[i, j]] == mode(kdirects[i], kdirects[j])
true

```

We also provide some shortcuts:

	kgridmap(kgrid::MonkhorstPack, mode = +)
	kgridmap(redkgrid::RedKgrid, mode = +)

"""
function kgridmap(kgrid::MonkhorstPack, mode = +)

	if !iszero(kgrid.kshift)
		error("Only work for Γ-centered kgrid!")
	end

	redkgrid = RedKgrid(kgrid)
	return kgridmap(redkgrid, redkgrid, mode)
end
function kgridmap(redkgrid::RedKgrid, mode = +)

	if !iszero(irredkgrid.kshift)
		error("Only work for Γ-centered kgrid!")
	end

	return kgridmap(redkgrid, redkgrid, mode)
end

function kgridmap(irredkgrid::IrredKgrid, mode = +)

	error("To be continued.")

	if !iszero(irredkgrid.kshift)
		error("Only work for Γ-centered kgrid!")
	end

	return kgridmap(irredkgrid, irredkgrid, mode)
end

function kgridmap(aimkdirects, kdirects, mode = +)

	Nk = length(kdirects)

	mapmatrix = Matrix{Int}(undef, Nk, Nk)
	Threads.@threads for I in CartesianIndices(mapmatrix)
		(i, j) = Tuple(I)
		Tk = mode(kdirects[i], kdirects[j])
		mapmatrix[I] = findfirst(k -> all(isinteger, k - Tk), aimkdirects)
	end

	return mapmatrix
end
