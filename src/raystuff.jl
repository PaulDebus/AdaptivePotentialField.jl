# SPDX-FileCopyrightText: 2023 Paul Debus <paul.debus@uni-weimar.de>
#
# SPDX-License-Identifier: GPL-3.0-or-later

function raycast_mesh_3d_cuda(mesh::Mesh, origin::Vector{Float32}, direction::Vector{Float32}, max_distance::Float32)
    vertices = CUDA.cu(mesh.vertices)
    faces = CUDA.cu(mesh.faces)
    origin = CUDA.cu(origin)
    direction = CUDA.cu(direction)
    intersects = CUDA.zeros(Float32, length(faces))
    
    CUDA.@allowscalar function raycast_mesh_3d_cuda_impl(vertices, faces, origin, direction, max_distance, intersects)
        index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        stride = gridDim().x * blockDim().x
        for i = index:stride:length(intersects)
            if i <= length(faces)
                face = faces[i]
                v1 = vertices[face[1]]
                v2 = vertices[face[2]]
                v3 = vertices[face[3]]
                intersects[i] = ray_triangle_intersect(origin, direction, v1, v2, v3, max_distance)
            end 
        end
        
    end
    
end # module AdaptivePotentialField


function ray_triangle_intersect(origin::Vector{Float32}, direction::Vector{Float32}, v1::Vector{Float32}, v2::Vector{Float32}, v3::Vector{Float32}, max_distance::Float32)
    e1 = v2 - v1
    e2 = v3 - v1
    pvec = cross(direction, e2)
    det = dot(e1, pvec)
    if det > -1e-8 && det < 1e-8
        return 0.0
    end
    inv_det = 1.0 / det
    tvec = origin - v1
    u = dot(tvec, pvec) * inv_det
    if u < 0.0 || u > 1.0
        return 0.0
    end
    qvec = cross(tvec, e1)
    v = dot(direction, qvec) * inv_det
    if v < 0.0 || u + v > 1.0
        return 0.0
    end
    t = dot(e2, qvec) * inv_det
    if t > 0.0 && t < max_distance
        return t
    end
    return 0.0
end
