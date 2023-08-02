# SPDX-FileCopyrightText: 2023 Paul Debus <paul.debus@uni-weimar.de>
#
# SPDX-License-Identifier: GPL-3.0-or-later

## imports
using AdaptivePotentialField
using Random
using Meshes
using Revise
using FileIO
using StaticArrays
using LinearAlgebra
using GeometryBasics: GeometryBasics
using Makie
using GLMakie
GLMakie.activate!()
Makie.inline!(false)


## load mesh and sample points
# geo_mesh = load("/home/paul/data/models/Scherkonde/pfeiler.ply")
geo_mesh = load("pfeiler.ply")
#geo_mesh = GeometryBasics.Cylinder3(GeometryBasics.Point(0,0,0), GeometryBasics.Point(0,0,20), 4.0)
mesh = convert(Meshes.SimpleMesh, geo_mesh)
sample_points, tris = random_points_and_normals(Random.GLOBAL_RNG, mesh, HomogeneousSampling(100000))
sample_points = [SVector{3,Float64}(coordinates(p)...) for p in sample_points];
normals = normal.(tris);
normals = [SVector{3, Float64}(n...) for n in normals];

## create a first viewpoint
coords = sample_points[1]
coords = coords + normals[1] * 4
view_dir = -normals[1]
vθ = Float64(asin(view_dir[3]))
vΦ = Float64(atan(view_dir[2], view_dir[1]))
v = PolarViewpoint(SVector{3, Float64}(coords), vθ, vΦ)
qv = QuatViewpoint(v)

##
weights = ones(length(sample_points))

@time potential(qv, sample_points, normals, weights)

## compute potential in a region around the viewpoint
function potential_region(vp::Viewpoint, points, normals, weights; d_max=2., steps=10, ang_max=deg2rad(20), ϵ_d=0.1, λ=0.5, ϵ=1e-2)
    if iseven(steps)
        @warn "grid size is even. Neutral movement is not contained"
    end
    xs = ys = zs = LinRange(-d_max, d_max, steps)
    θs = Φs = LinRange(-ang_max, ang_max, steps)
    pot = zeros(length(xs), length(ys), length(zs), length(θs), length(Φs))
    for (i, x) in enumerate(xs)
        for (j, y) in enumerate(ys)
            for (k, z) in enumerate(zs)
                for (l, θ) in enumerate(θs)
                    for (m, Φ) in enumerate(Φs)
                        mov = PolarMovement(SVector{3, Float64}([x, y, z]), θ, Φ)
                        v = apply(vp, mov)
                        pot[i,j,k,l,m] = potential(v, points, normals, weights, ϵ_d=ϵ_d)
                    end
                end
            end
        end
    end
    pot
end
function scale_potential_horizontal!(pot::Array, λ=0.7)
    # to prefer horizontal movements, scale movement along z axis by something <1
    center = Int(ceil(size(pot,3)/2))
    θ0 = ϕ0 = Int(ceil(size(pot,4)/2))
    for idx in eachindex(IndexCartesian(), pot)
        if idx[3] != center || idx[4] != θ0
            pot[idx] *= λ 
        elseif idx[5] != ϕ0 && (abs(idx[2] - center) + abs(idx[1] - center)) > 2
            pot[idx] *= 1
        end
        # pot[idx] *= 1. - abs(idx[3] - center) / 10.
    end
end

function index_to_movement(idx::CartesianIndex; d_max=2., steps=10, ang_max=deg2rad(20))
    xs = ys = zs = LinRange(-d_max, d_max, steps)
    θs = Φs = LinRange(-ang_max, ang_max, steps)
    x = xs[idx[1]]
    y = ys[idx[2]]
    z = zs[idx[3]]
    θ = θs[idx[4]]
    Φ = Φs[idx[5]]
    PolarMovement(SVector{3, Float64}([x, y, z]), θ, Φ)
end
##
@time pot = potential_region(qv, sample_points, normals, weights, steps=5);
scale_potential_horizontal!(pot)
## 
argmax(pot)
index_to_movement(CartesianIndex(3, 5, 3, 3, 3), steps=5)

## try planning
vps = [qv]
moves = []
steps = 5
d_max = 1.5
ϵ_d = 0.01
weights = ones(length(sample_points));
##
@time for i in 1:50
    update_weights!(weights, sample_points, normals, vps[end], ϵ_d=ϵ_d)
    region = potential_region(vps[end], sample_points, normals, weights, 
                                steps=steps, d_max=d_max, ϵ_d=ϵ_d)
    scale_potential_horizontal!(region, 0.4)
    highest_idx = argmax(region)
    move = index_to_movement(highest_idx, steps=steps, d_max=d_max)
    push!(moves, move)
    push!(vps, apply(vps[end], move))
    # visualization
    delete!(scene, viewdir_vis)
    vp_pos = [GeometryBasics.Point3(v.position...) for v in vps]
    viewdirs = GeometryBasics.Vec3.(viewdir.(vps))
    viewdir_vis = arrows!(scene, vp_pos, viewdirs.*4, color=:red)
    delete!(scene, point_vis)
    point_vis = scatter!(scene, sample_points, color=weights)
    
end

##
scene = Scene()
cam3d!(scene)
Makie.mesh!(scene, geo_mesh)
center!(scene)
##
# delete!(scene, viewdir_vis)
vp_pos = [GeometryBasics.Point3(v.position...) for v in vps]
viewdirs = GeometryBasics.Vec3.(viewdir.(vps))
viewdir_vis = arrows!(scene, vp_pos, viewdirs.*4, color=:red)
# delete!(scene, point_vis)
point_vis = scatter!(scene, sample_points, color=weights) 
#delete!(scene, vps_vis)
# vps_vis = meshscatter!(scene, vp_pos, markersize=.2)
# delete!(scene, path_vis)
# vp_pos = [GeometryBasics.Point3(v.position...) for v in vps]
# path_vis = lines!(scene, vp_pos, color=:red)

###########
# TODO: continue here
# - check angle constraint for visibility
# - use ImageProjectiveGeometry.Camera instead of Viewpoint? Could make things easier
# - integrate all parameters and data into one "scene" object (need better name)
# - constrain the movement to more consistent (depends on if the rest already works)
# - use bppc for weight update


## find out how to convert from viewpoint space to global space
function to_global(v::Viewpoint, mov::QuatMovement)
    trans_global = AdaptivePotentialField.rotate(mov.direction, v.orientation)
    rot_global = mov.rotation
    QuatMovement(SVector(trans_global...), rot_global)
end
# qq = QuatViewpoint(PolarViewpoint(SVector{3, Float64}([0., 0., 1.]), 0., pi/4))
qv2 = QuatViewpoint(PolarViewpoint(qv.position, -pi/4, pi/1))
qm = QuatMovement(PolarMovement(SVector{3, Float64}([0., -1., 0.]), pi/8., -0pi/5))
gm =to_global(qv, qm)
pp = PolarViewpoint(apply(qv2, qm))
println(qv.position)
println(pp.position)

#qv2 = PolarViewpoint(qv.position, -pi/4, pi/1)

delete!(scene, viewdir_vis)
vps = [qv2, pp]
vp_pos = [GeometryBasics.Point3(v.position...) for v in vps]
viewdirs = GeometryBasics.Vec3.(viewdir.(vps))
viewdir_vis = arrows!(scene, vp_pos, viewdirs.*4, color=[:red, :green])