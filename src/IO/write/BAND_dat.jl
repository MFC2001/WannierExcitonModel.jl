function Base.write(io::IO, band::BSE_uEigen, ::Type{BAND_dat}; kwargs...)
    bandmatrix = _eigen2val(band.eigen)
    write(io, bandmatrix, BAND_dat; kwargs...)
end
function Base.write(io::IO, band::AbstractVector{<:BSE_uEigen}, ::Type{BAND_dat}; kwargs...)
    bandmatrix = _eigen2val(map(x->x.eigen, band))
    write(io, bandmatrix, BAND_dat; kwargs...)
end
