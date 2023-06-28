module AdaptivePotentialField

using CUDA
using Meshes
using Random
using GeometryBasics: GeometryBasics
using FileIO
using ImageProjectiveGeometry
using Quaternionic
using LinearAlgebra, StaticArrays



include("utils.jl")
include("viewpoints.jl")
include("visualization.jl")

function update_weights!(weights::Vector, vp::Viewpoint; ϵ_d=0.1, λ=0.5, ϵ=1e-2)
    for i in 1:length(sample_points)
        if is_point_visible_float(sample_points[i], normals[i], vp.position, viewdir(vp), 4., π/3, cos(π/4), ϵ_d) > ϵ
            weights[i] *= λ
        end
    end
end

function potential(x::Viewpoint; target_distance=4., fov=π/3, ϵ_d=0.1, angle_tol = π/4, make_training_data=false)
    cos_max = cos(angle_tol)
    pos = x.position
    dir = viewdir(x)
    
    ThreadsX.mapreduce(+, sample_points, normals, weights) do point, normal, weight
        vis = is_point_visible_float(point, normal, pos, dir, target_distance, fov, cos_max, ϵ_d)
        return vis *  weight
    end
end


end