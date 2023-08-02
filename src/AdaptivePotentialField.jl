# SPDX-FileCopyrightText: 2023 Paul Debus <paul.debus@uni-weimar.de>
#
# SPDX-License-Identifier: GPL-3.0-or-later

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
include("surface_points.jl")
export update_weights!, potential, is_point_visible_float

# soft isless function for smooth transitions and better gradients
@inline cont_isless(x, y; a=100) = 1 / (1+exp(-(y-x)*a))

function is_point_visible_float(point::SurfacePoint, position, viewdir, target_distance, fov=π/3, 
                                cos_min=cos(π/4), ϵ_d=0.1)
    # approximate the computation of the visibility of a surface point from a viewpoint
    # based on the distance to the point, the field of view and the angle between the normal 
    # and the view direction. does not check for true visibility with collision checks.
    
    # admissible distance range: is the distance between the point and the viewpoint within
    # the target distance +- ϵ_d defined by the resolution requirements?
    dist = norm(point.point - position)
    δ = abs(dist - target_distance) / target_distance # relative distance to target distance
    d_vis = cont_isless(δ, ϵ_d, a=1000) # use stronger exponent since we have small values
    
    # check if the surface point is within the field of view by checking the angle between
    # the view direction and the vector from the viewpoint to the surface point (viewray)
    # using spherical view, not a rectangular image
    viewray = (point.point - position)./dist
    d_fov = cont_isless(cos(fov/2), viewdir ⋅ viewray)
    
    # check the incidence angle of the viewing ray on the surface to avoid 
    # shallow angles that lead to strong distortions
    d_ang = cont_isless(cos_min, viewray ⋅ -point.normal)

    # return the product of the three factors to conserve zeros
    return d_vis * d_fov * d_ang
end

function update_weights!(points::Vector{SurfacePoint}, vp::Viewpoint; target_distance=4, fov=π/3, ϵ_d=0.1, angle_tol = π/4, λ=0.5, ϵ=1e-2)
    # update the weights of the surface points based if they are visibile from the viewpoint
    for point in points
        if is_point_visible_float(point, vp.position, viewdir(vp), target_distance, fov, cos(angle_tol), ϵ_d) > ϵ
            reduce_weight!(point, λ)
        end
    end
end

function potential(x::Viewpoint, points::Vector{SurfacePoint}; target_distance=4., fov=π/3, ϵ_d=0.1, angle_tol = π/4, safety_distance=1.) where T
    # the potential of a viewpoint candidate is the weighted sum of the visibilities
    # of all surface points by the weight of each point
    cos_max = cos(angle_tol)
    pos = x.position
    dir = viewdir(x)

    if any(distance.(points, pos) .< safety_distance)
        return zero(T)
    end
    
    # TODO: is ThreadsX faster than normal vectorization?
    ThreadsX.mapreduce(+, points) do point
        vis = is_point_visible_float(point, pos, dir, target_distance, fov, cos_max, ϵ_d)
        return vis *  point.weight
    end
end


end