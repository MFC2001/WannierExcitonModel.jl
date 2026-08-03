
function ChernNumber(::Union{Val{:QuantumGeometry}, Val{:QG}}, kgrid, TB::AbstractTightBindModel; bandindex = 1)

	BC = BerryCurvature(Val(:QuantumGeometry), kgrid, TB, bandindex)

	𝐚, 𝐛, 𝐜 = basisvectors(reciprocal(TB.lattice))

	if length(bandindex) == 1
		if kgrid.kgrid_size[1] == 1
			S = abs((𝐛×𝐜)[1])
			BC = BC[2, 3, :]
		elseif kgrid.kgrid_size[2] == 1
			S = abs((𝐚×𝐜)[2])
			BC = BC[3, 1, :]
		elseif kgrid.kgrid_size[3] == 1
			S = abs((𝐚×𝐛)[3])
			BC = BC[1, 2, :]
		else
			error("Only calculate Chern Number for 2D.")
		end
		return sum(BC) * S / length(kgrid) / 2π
	else
		if kgrid.kgrid_size[1] == 1
			S = abs((𝐛×𝐜)[1])
			BC = BC[:, :, 2, 3, :]
		elseif kgrid.kgrid_size[2] == 1
			S = abs((𝐚×𝐜)[2])
			BC = BC[:, :, 3, 1, :]
		elseif kgrid.kgrid_size[3] == 1
			S = abs((𝐚×𝐛)[3])
			BC = BC[:, :, 1, 2, :]
		else
			error("Only calculate Chern Number for 2D.")
		end
		TBC = similar(BC, size(BC, 3))
		Threads.@threads for k in axes(BC, 3)
			TBC[k] = tr(BC[:, :, k])
		end
		return sum(TBC) * S / length(kgrid) / 2π
	end
end
function ChernNumber(::Union{Val{:WilsonLoop}, Val{:WL}}, kgrid::MonkhorstPack, TB::AbstractTightBindModel; bandindex = 1)
	(center, BC) = BerryCurvature(Val(:WilsonLoop), kgrid, TB; bandindex)
	return ChernNumber(BC, kgrid, TB.lattice)
end
