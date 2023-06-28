using StaticArrays
using Makie
using GLMakie
using ImageProjectiveGeometry

using AdaptivePotentialField

vp = PolarViewpoint(SVector(0., 0, 0), 0., 0.)
scene = Scene()
cam3d!(scene)

cam = convert(ImageProjectiveGeometry.Camera, vp)

AdaptivePotentialField.plotcamera(cam, 4., scene)