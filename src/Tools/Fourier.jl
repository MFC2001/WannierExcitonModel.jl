function Fourier(fk::AbstractVector{<:Number}, kmesh::KData; coe=kmesh.N, mode=+)
    Nmesh = kmesh.Nmesh
    meshindex = kmesh.meshindex

    fR = Vector{ComplexF64}(undef, kmesh.N)
    Threads.@threads for i in eachindex(fR)
        R = mode.(0, meshindex[i] .// Nmesh)
        fR[i] = sum(fk .* map(k -> cispi(2 * k ⋅ R), meshindex))
    end

    fR ./= coe

    return fR
end