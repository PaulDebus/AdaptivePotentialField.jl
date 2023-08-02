# SPDX-FileCopyrightText: 2023 Paul Debus <paul.debus@uni-weimar.de>
#
# SPDX-License-Identifier: GPL-3.0-or-later

## 
using AdaptivePotentialField
using Plots 
import Plots: plot, plot!
using LinearAlgebra
using StaticArrays

## some helper classes
struct Box{T}
    x::T
    y::T
    w::T
    h::T
end
struct Point{T}
    x::T
    y::T
end
function plot(b::Box)
    coords = [(b.x, b.y), (b.x+b.w, b.y), (b.x+b.w, b.y+b.h), (b.x, b.y+b.h), (b.x,b.y)]
    plot(Shape(coords), label="", c=:black)
end
function plot!(b::Box)
    coords = [(b.x, b.y), (b.x+b.w, b.y), (b.x+b.w, b.y+b.h), (b.x, b.y+b.h), (b.x,b.y)]
    plot!(Shape(coords), label="", c=:black)
end
function dist(b::Box, point::Point)
    cx = b.x + 0.5b.w
    cy = b.y + 0.5b.h
    dx = abs(cx - point.x) - (b.w * 0.5)
    dy = abs(cy - point.y) - (b.h * 0.5)
    return sqrt((dx * (dx > 0))^2 + (dy * (dy > 0))^2)
end

## 
b1 = Box(10,10, 70,10)
b2 = Box(80,10, 10, 50)
boxes = [b1, b2]


##
x = LinRange(41,46, 200)
y = LinRange(-3,12, 300)
points = [SVector(45., i, 0.) for i in LinRange(1, 5, 10)]
normals = [SVector(-1., 0., 0.) for i in LinRange(1, 5, 10)]
weights = ones(length(points))
##
weights[1:6] .*= 0.5
c = [potential(from_vector(QuatViewpoint, [X, Y, 0., 0., 0.]), points, normals, weights, target_distance=2., ϵ_d=0.2) for X in x, Y in y]

##
colors = cgrad([:white, :white, :green], [0., 0.25, 1.])

heatmap(x, y, c', c=colors,  dpi=300)
weight_colors = cgrad([:white, :black], [0., 1.])
scatter!([p[1] for p in points], [p[2] for p in points], zcolor=weights.*5., label="", m=weight_colors, msize=8)

##
savefig("/home/paul/work/cloud/data_sync/research/publications/2023_ISPRS_GeospatialWeek/figures/update2_2.png")