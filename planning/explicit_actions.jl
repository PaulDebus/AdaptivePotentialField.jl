# SPDX-FileCopyrightText: 2023 Paul Debus <paul.debus@uni-weimar.de>
#
# SPDX-License-Identifier: GPL-3.0-or-later

##
using AdaptivePotentialField
using Random
using Meshes
using Revise
using FileIO
using StaticArrays
using LinearAlgebra
using GeometryBasics: GeometryBasics
using Statistics
using Quaternionic
using Makie
using GLMakie
GLMakie.activate!()
Makie.inline!(false)

##
geo_mesh = load("/home/paul/work/cloud/data_sync/research/bppc/Artifical-House.ply")
# geo_mesh = load("pfeiler.ply")
mesh = convert(Meshes.SimpleMesh, geo_mesh)
sample_points, tris = random_points_and_normals(Random.GLOBAL_RNG, mesh, HomogeneousSampling(100000))
sample_points = [SVector{3,Float64}(coordinates(p)...) for p in sample_points];
normals = normal.(tris);
normals = [SVector{3, Float64}(n...) for n in normals];
ground_level = minimum(getindex.(sample_points, 3))
# sample_points = [p for p in sample_points if p[3] > ground_level + 0.01]
# normals = [n for (n, p) in zip(normals, sample_points) if p[3] > ground_level + 0.01]

## create a first viewpoint
function random_viewpoint(points, point_normals, dist=4.)
    idx = rand(eachindex(points))
    println("pos: ", points[idx])
    println("norm: ", point_normals[idx])
    view_dir = -point_normals[idx]
    view_dir = view_dir / norm(view_dir)
    coords = points[idx] + (-view_dir * dist)
    vθ = Float64(asin(view_dir[3]))
    vΦ = Float64(atan(view_dir[2], view_dir[1]))
    v = PolarViewpoint(SVector{3, Float64}(coords), vθ, vΦ)
    QuatViewpoint(v), idx
end
vp = random_viewpoint(sample_points, normals)

##
d_target = 4.
fov = deg2rad(60)
d_max = 2*tan(fov/2)*d_target / 3 ## overlap 66%
rot_max = deg2rad(20)
ϵ_d = 0.1
rc = cos(rot_max)
rs = sin(rot_max)
admissible_local_movements = QuatMovement.([
    PolarMovement(SVector{3, Float64}([0, 0, d_max]), 0., 0.),
    PolarMovement(SVector{3, Float64}([0, 0, -d_max]), 0., 0.),
    PolarMovement(SVector{3, Float64}([0, d_max, 0]), 0., 0.),
    PolarMovement(SVector{3, Float64}([0, -d_max, 0]), 0., 0.),
    PolarMovement(SVector{3, Float64}([d_max, 0, 0]), 0., 0.),
    PolarMovement(SVector{3, Float64}([-d_max, 0, 0]), 0., 0.),
    PolarMovement(SVector{3, Float64}([d_target*(1-rc), 0., d_target*rs]), rot_max, 0.), # rotate up
    PolarMovement(SVector{3, Float64}([d_target*(1-rc), 0., -d_target*rs]), -rot_max, 0.), # rotate down
    PolarMovement(SVector{3, Float64}([d_target*(1-rc), d_target*rs, 0]), 0., -rot_max), # rotate left
    PolarMovement(SVector{3, Float64}([d_target*(1-rc), -d_target*rs, 0]), 0., rot_max), # rotate right
])
admissible_local_movements_polar = [
    PolarMovement(SVector{3, Float64}([0, 0, d_max]), 0., 0.),
    PolarMovement(SVector{3, Float64}([0, 0, -d_max]), 0., 0.),
    PolarMovement(SVector{3, Float64}([0, d_max, 0]), 0., 0.),
    PolarMovement(SVector{3, Float64}([0, -d_max, 0]), 0., 0.),
    PolarMovement(SVector{3, Float64}([d_max, 0, 0]), 0., 0.),
    PolarMovement(SVector{3, Float64}([-d_max, 0, 0]), 0., 0.),
    PolarMovement(SVector{3, Float64}([d_target*(1-rc), 0., d_target*rs]), rot_max, 0.), # rotate up
    PolarMovement(SVector{3, Float64}([d_target*(1-rc), 0., -d_target*rs]), -rot_max, 0.), # rotate down
    PolarMovement(SVector{3, Float64}([d_target*(1-rc), d_target*rs, 0]), 0., -rot_max), # rotate left
    PolarMovement(SVector{3, Float64}([d_target*(1-rc), -d_target*rs, 0]), 0., rot_max), # rotate right
]
move_weights = [0.5, 0.5, 1., 1., 0., 0., 0.5, 0.5, 1., 1.] # prefer horizontal movements
function to_global(v::Viewpoint, mov::QuatMovement)
    turn_back = 
    trans_global = AdaptivePotentialField.rotate(mov.direction, v.orientation)

    rot_global = mov.rotation # TODO: rotate around global z, right now is rotated z
    QuatMovement(SVector(trans_global...), rot_global)
end

##
function global_apply(vp::QuatViewpoint, mov::PolarMovement) 
    trans_global = SVector(AdaptivePotentialField.rotate(mov.direction, vp.orientation)...)
    QuatViewpoint(trans_global+vp.position, exp(imz*mov.dΦ/2) *vp.orientation * exp(imy*mov.dΘ/2))
end

##
scene = Scene()
cam3d!(scene)
Makie.mesh!(scene, geo_mesh)
center!(scene)

##
qv, idx = random_viewpoint(sample_points, normals, d_target)
weights[idx] = 0.
vps = [qv]
vp_pos = [GeometryBasics.Point3(v.position...) for v in vps]
viewdirs = GeometryBasics.Vec3.(viewdir.(vps))
delete!(scene, viewdir_vis)
delete!(scene, point_vis)
viewdir_vis = arrows!(scene, vp_pos, viewdirs.*4, color=:red)
point_vis = scatter!(scene, sample_points, color=weights) 


##

flightpath = [qv]
weights = ones(length(sample_points))
meanweights = [1.]
##
for i in 1:250
    update_weights!(weights, sample_points, normals, flightpath[end], ϵ_d=ϵ_d)
    candidates =  [global_apply(flightpath[end], mov) for mov in admissible_local_movements_polar]
    potentials = [potential(v, sample_points, normals, weights) for v in candidates]
    potentials .*= move_weights.^1.5
    push!(flightpath, candidates[argmax(potentials)])

    vp_pos = [GeometryBasics.Point3(v.position...) for v in flightpath]
    viewdirs = GeometryBasics.Vec3.(viewdir.(flightpath))
    delete!(scene, viewdir_vis)
    delete!(scene, point_vis)
    viewdir_vis = arrows!(scene, vp_pos, viewdirs.*4, color=:red)
    point_vis = scatter!(scene, sample_points, color=weights) 
    if i % 10 == 0
        push!(meanweights, mean(weights))
    end
end

##
for i in 1:50
    update_weights!(weights, sample_points, normals, flightpath[end], ϵ_d=ϵ_d)
    admissible_movements = [to_global(flightpath[end], mov) for mov in admissible_local_movements]
    candidates =  [apply(flightpath[end], mov) for mov in admissible_movements]
    potentials = [potential(v, sample_points, normals, weights) for v in candidates]
    potentials .*= move_weights.^1.5
    push!(flightpath, candidates[argmax(potentials)])

    vp_pos = [GeometryBasics.Point3(v.position...) for v in flightpath]
    viewdirs = GeometryBasics.Vec3.(viewdir.(flightpath))
    delete!(scene, viewdir_vis)
    delete!(scene, point_vis)
    viewdir_vis = arrows!(scene, vp_pos, viewdirs.*4, color=:red)
    point_vis = scatter!(scene, sample_points, color=weights) 
    if i % 10 == 0
        push!(meanweights, mean(weights))
    end
end

