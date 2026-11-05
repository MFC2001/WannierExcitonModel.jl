function Base.write(io::IO, band::BSE_uEigen, ::Type{wannier90_eig}; kwargs...)
    bandmatrix = _eigen2val(band.eigen)
    write(io, bandmatrix, wannier90_eig; kwargs...)
end
function Base.write(io::IO, band::AbstractVector{<:BSE_uEigen}, ::Type{wannier90_eig}; kwargs...)
    bandmatrix = _eigen2val(map(x->x.eigen, band))
    write(io, bandmatrix, wannier90_eig; kwargs...)
end
