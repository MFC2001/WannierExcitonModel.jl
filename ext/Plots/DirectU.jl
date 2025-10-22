function Plots.plot(U::HR, TB::AbstractTightBindModel; kwargs...)
	return Plots.plot(U, TB.lattice, TB.orb_location; kwargs...)
end
function Plots.plot(U::HR, lattice::Lattice, orblocat::AbstractVector{<:AbstractVector}; ϵ::Real = 1,
	_value = real.(U.value),
	ylims = (-0.1, maximum(_value) + 0.5))

	r_point = Vector{Float64}(undef, length(_value))
	for i in eachindex(_value)
		r_point[i] = norm(lattice * (U.path[i, 1:3] + orblocat[U.path[i, 5]] - orblocat[U.path[i, 4]]))
	end

	r = range(0.01, maximum(r_point) + 5, length = 101)
	V = RealInverseR(; ϵ)
	V = V.(r)

	p = Plots.plot(r, V;
		title = "",
		linewidth = 1.2,
		linecolor = :black,
		xlims = (0, r[end]),
		ylims,
		legend = false,
		size = (800, 600),
	)
	Plots.plot!(r_point, _value;
		seriestype = :scatter,
	)

	return p
end
