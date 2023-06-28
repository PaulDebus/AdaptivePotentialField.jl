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
geo_mesh = load("pfeiler.ply")
geo_mesh = GeometryBasics.Cylinder3(GeometryBasics.Point(0,0,0), GeometryBasics.Point(0,0,20), 4.0)
mesh = convert(Meshes.SimpleMesh, geo_mesh)
sample_points, tris = random_points_and_normals(Random.GLOBAL_RNG, mesh, HomogeneousSampling(100000))
sample_points = [SVector{3,Float64}(coordinates(p)...) for p in sample_points];
normals = normal.(tris);
normals = [SVector{3, Float64}(n...) for n in normals];

## create a first viewpoint
coords = sample_points[1]
coords = coords + normals[1] * 4
vΦ = Float64(asin(coords[3] / norm(coords)))
vθ = Float64(atan(coords[2], coords[1]))
v = PolarViewpoint(SVector{3, Float64}(coords), vθ, vΦ)
qv = QuatViewpoint(v)

##
weights = ones(length(sample_points))

@time potential(qv, sample_points, normals, weights)

## compute potential in a region around the viewpoint
function potential_region(vp::Viewpoint, points, normals, weights; d_max=2., steps=10, ang_max=deg2rad(20))
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
                        pot[i,j,k,l,m] = potential(v, points, normals, weights)
                    end
                end
            end
        end
    end
    pot
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
@time pot = potential_region(qv, sample_points, normals, weights, steps=7);

## 
index_to_movement(argmax(pot), steps=7)

##
function update_weights!(weights::Vector, sample_points, normals, vp::Viewpoint; ϵ_d=0.1, λ=0.5, ϵ=1e-2)
    for i in eachindex(sample_points)
        if is_point_visible_float(sample_points[i], normals[i], vp.position, viewdir(vp), 4., π/3, cos(π/4), ϵ_d) > ϵ
            weights[i] *= λ
        end
    end
end

## try planning
vps = [qv]
moves = []
steps = 5
d_max = 1.5
weights = ones(length(sample_points));
##
@time for i in 1:25
    update_weights!(weights, sample_points, normals, vps[end])
    region = potential_region(vps[end], sample_points, normals, weights, steps=steps, d_max=d_max)
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
delete!(scene, viewdir_vis)
vp_pos = [GeometryBasics.Point3(v.position...) for v in vps]
viewdirs = GeometryBasics.Vec3.(viewdir.(vps))
viewdir_vis = arrows!(scene, vp_pos, viewdirs.*4, color=:red)
delete!(scene, point_vis)
point_vis = scatter!(scene, sample_points, color=weights) 
#delete!(scene, vps_vis)
# vps_vis = meshscatter!(scene, vp_pos, markersize=.2)

###########
# TODO: continue here
# - make fov smaller
# - check angle constraint for visibility
# - why are the viewpoints moving closer to the surface???
# - integrate all parameters and data into one "scene" object (need better name)
# - constrain the movement to more consistent (depends on if the rest already works)
# - use bppc for weight update