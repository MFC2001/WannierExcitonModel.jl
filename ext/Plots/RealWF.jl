@userplot ColorScatter

@recipe function f(cs::ColorScatter; normal = nothing)
    values, locations = cs.args
    pts = map(r -> Vec3{Float64}(r), locations)

    seriestype := :scatter
    marker_z   := values

    if isnothing(normal)
        # 三维分支
        xs = getindex.(pts, 1)
        ys = getindex.(pts, 2)
        zs = getindex.(pts, 3)
        return (xs, ys, zs)
    else
        # 投影分支
        frame = plane_frame(normal)
        pu = map(r -> dot(r, frame.e1), pts)
        pv = map(r -> dot(r, frame.e2), pts)
        return (pu, pv)
    end
end


# 色度散点图
function colorscatter(values::AbstractVector{<:Real},
    locations::AbstractVector{<:AbstractVector{<:Real}};
    normal=nothing,
    kwargs...,
)
    locations = map(r->Vec3{Float64}(r), locations)
    if isnothing(normal)
        x = map(r->r[1], locations)
        y = map(r->r[2], locations)
        z = map(r->r[3], locations)

        #三维散点图
    else
        frame = plane_frame(normal)
        pu = map(locations) do r
            r ⋅ frame.ê₁
        end
        pv = map(locations) do r
            r ⋅ frame.ê₂
        end

        # 二维散点图
    end
    return p
end
