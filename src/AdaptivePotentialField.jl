module AdaptivePotentialField

using Meshes
using Random
using GeometryBasics: GeometryBasics
using FileIO
using ImageProjectiveGeometry
using Quaternionic
using LinearAlgebra, StaticArrays
using ThreadsX



include("utils.jl")
include("viewpoints.jl")
include("visualization.jl")
export update_weights!, potential, is_point_visible_float

@inline cont_isless(x, y; a=100) = 1 / (1+exp(-(y-x)*a))
function is_point_visible_float(point, normal, position, viewdir, target_distance, fov=π/3, cos_min=cos(π/4), ϵ_d=0.1)
    dist = norm(point - position)
    δ = abs(dist - target_distance) / target_distance 
    
    # d_vis = 1. - clamp((δ / ϵ_d)^10, 0, 1) # admissible distance range
    d_vis = cont_isless(δ, ϵ_d/2) # TODO: is this a bad hack?
    # d_vis = δ < ϵ_d
    
    viewray = (point - position)./dist
    #d_fov = clamp(((viewray ⋅ viewdir) / cos(fov)), 0., 1.)^10 # field of view
    d_fov = cont_isless(cos(fov/2), viewdir ⋅ viewray)
    
    #d_ang = clamp(((viewray ⋅ -normal) / cos_min), 0., 1.)^10 # angle between normal and viewdir
    d_ang = cont_isless(cos_min, viewray ⋅ -normal)
    return d_vis * d_fov * d_ang
end

function update_weights!(weights::Vector, points, normals, vp::Viewpoint; target_distance=4, fov=π/3, ϵ_d=0.1, angle_tol = π/4, λ=0.5, ϵ=1e-2)
    for i in eachindex(points)
        if is_point_visible_float(points[i], normals[i], vp.position, viewdir(vp), target_distance, fov, cos(angle_tol), ϵ_d) > ϵ
            weights[i] *= λ
        end
    end
end

function potential(x::Viewpoint, points::Vector{SVector{3,T}}, normals::Vector{SVector{3,T}}, weights::Vector; target_distance=4., fov=π/3, ϵ_d=0.1, angle_tol = π/4) where T
    cos_max = cos(angle_tol)
    pos = x.position
    dir = viewdir(x)

    safety_distance = 1.
    if any([norm(p - pos) for p in points] .< safety_distance)
        return zero(T)
    end
    
    ThreadsX.mapreduce(+, points, normals, weights) do point, normal, weight
        vis = is_point_visible_float(point, normal, pos, dir, target_distance, fov, cos_max, ϵ_d)
        return vis *  weight
    end
end


end