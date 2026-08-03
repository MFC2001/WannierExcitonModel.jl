
function BAND_hrh(kpoints, hrh, ::Val{true})
	band = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, length(kpoints))
	Threads.@threads for i in eachindex(kpoints)
		band[i] = eigen!(hrh(kpoints[i]))
	end
	return band
end
function BAND_hrh(kpoints, hrh, ::Val{false})
	band = Matrix{Float64}(undef, numorb(hrh), length(kpoints))
	Threads.@threads for i in eachindex(kpoints)
		band[:, i] = eigvals!(hrh(kpoints[i]))
	end
	return band
end
function BAND_hrh(kpoints, hrh, orblocat, ::Val{true})
	band = Vector{Eigen{ComplexF64, Float64, Matrix{ComplexF64}, Vector{Float64}}}(undef, length(kpoints))
	Threads.@threads for i in eachindex(kpoints)
		band[i] = eigen!(hrh(kpoints[i], orblocat))
	end
	return band
end
function BAND_hrh(kpoints, hrh, orblocat, ::Val{false})
	return BAND_hrh(kpoints, hrh, Val(false))
end
