# SPDX-FileCopyrightText: 2023 Paul Debus <paul.debus@uni-weimar.de>
#
# SPDX-License-Identifier: GPL-3.0-or-later

using StaticArrays
using LinearAlgebra

export SurfacePoint, distance, reduce_weight!

mutable struct SurfacePoint{T<:Real}
    const point::SVector{3, T}
    const normal::SVector{3, T}
    weight::Real
end
function SurfacePoint(point, normal; init=1.)
    SurfacePoint(point, normal, init)
end

distance(p::SurfacePoint, q::SVector) = norm(p.point - q)
distance(p::SurfacePoint, vp::Viewpoint) = distance(p, vp.position)

function reduce_weight!(point::SurfacePoint, λ=0.5)
    point.weight *= λ
end