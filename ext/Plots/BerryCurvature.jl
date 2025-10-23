
function Plots.plot(BC::BerryCurvature_wilsonloop;
		layout = Plots.@layout [a{0.95w} b])
	kgrid_size = BC.center.kgrid_size
	if kgrid_size[1] == 1
		kx = 2
		ky = 3
	elseif kgrid_size[2] == 1
		kx = 3
		ky = 1
	elseif kgrid_size[3] == 1
		kx = 1
		ky = 2
	end
	rlattice = reciprocal(BC.lattice)
	center_car = map(k -> rlattice * k, BC.center)
	kx = map(k -> k[kx], center_car)
	ky = map(k -> k[ky], center_car)

	if BC.buffer isa AbstractVector
		BC_v = BC.buffer * BC.δS / 2π
	elseif BC.buffer isa AbstractArray{<:Number, 3}
		BC_v = Vector{Float64}(undef, size(BC.buffer, 3))
		Threads.@threads for k in axes(BC.buffer, 3)
			BC_v[k] = real(tr(view(BC.buffer, :, :, k))) * BC.δS / 2π
		end
	else
		error("Wrong BerryCurvature!")
	end

	x_range = extrema(kx)
	y_range = extrema(ky)
	common_range = (min(x_range[1], y_range[1]), max(x_range[2], y_range[2]))
	δ = (common_range[2] - common_range[1]) / 20
	common_range = (common_range[1] - δ, common_range[2] + δ)

	p = Plots.plot(kx, ky;
		size = (600, 600),
		aspect_ratio = :equal,
		xlims = common_range,
		ylims = common_range,
		legend = false,
		seriestype = :scatter,
		marker_z = BC_v,
		markersize = 3,
		markerstrokewidth = 0,
		color = :viridis,
		colorbar = false,
		xlabel = "kx",
		ylabel = "ky",
		title = "")
	return p
end
