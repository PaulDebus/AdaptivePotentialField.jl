# SPDX-FileCopyrightText: 2023 Paul Debus <paul.debus@uni-weimar.de>
#
# SPDX-License-Identifier: GPL-3.0-or-later

using StaticArrays
using Makie
using GLMakie
using ImageProjectiveGeometry
GLMakie.activate!()
Makie.inline!(false)
using AdaptivePotentialField
using GeometryBasics

vp = PolarViewpoint(SVector(0., 0, 0), 0., 0.)
qp = QuatViewpoint(vp)
scene = Scene()
cam3d!(scene)
arrows!(scene, [Point3(0, 1., 0)], [Vec3(1., 0, 0)], color=:red)
scene

cam = convert(ImageProjectiveGeometry.Camera, qp)
cam.rows = 100
cam.cols = 100

subscene = AdaptivePotentialField.plotcamera(cam, 4., scene);
center!(scene)