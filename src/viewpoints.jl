# SPDX-FileCopyrightText: 2023 Paul Debus <paul.debus@uni-weimar.de>
#
# SPDX-License-Identifier: GPL-3.0-or-later

export Viewpoint, PolarViewpoint, QuatViewpoint
export viewdir, to_vector, from_vector
export Movement, PolarMovement, QuatMovement
export apply, isvalid, make_valid, interpolate


# utils
function viewdir_from_spherical_x_coordinates(θ, ϕ)
    return [cos(ϕ)*cos(θ), sin(ϕ)*cos(θ),  -sin(θ)]
end
rotate(p, q::Rotor) = (q*QuatVec(p)*conj(q)).im


############################################################
# Viewpoint
############################################################
abstract type Viewpoint end
struct PolarViewpoint{T} <: Viewpoint
    position::SVector{3, T}
    θ::T # pitch, look up/down
    Φ::T # zaw, rotate around up==z
end
struct QuatViewpoint{T} <: Viewpoint
    position::SVector{3, T}
    orientation::Rotor
end

# get normalized view direction
viewdir(vp::PolarViewpoint) = viewdir_from_spherical_x_coordinates(vp.θ, vp.Φ)
viewdir(vp::QuatViewpoint) = rotate([1., 0., 0.], vp.orientation)

# convert between viewpoint types
QuatViewpoint(pv::PolarViewpoint) = QuatViewpoint(pv.position, from_euler_angles(pv.Φ,pv.θ,0))
function PolarViewpoint(qv::QuatViewpoint)
    eulers = to_euler_angles(qv.orientation)
    PolarViewpoint(qv.position, eulers[2], eulers[1] + eulers[3])
end
Base.convert(::Type{QuatViewpoint}, pv::PolarViewpoint) = QuatViewpoint(pv)
Base.convert(::Type{PolarViewpoint}, qv::QuatViewpoint) = PolarViewpoint(qv)

# convert to/from vectors for serialization and interaction with autodiff
to_vector(qv::QuatViewpoint) = to_vector(PolarViewpoint(qv))
to_vector(pv::PolarViewpoint) = [pv.position..., pv.θ, pv.Φ]
from_vector(::Type{QuatViewpoint}, v) = QuatViewpoint(from_vector(PolarViewpoint, v))
from_vector(::Type{PolarViewpoint}, v) = PolarViewpoint(SVector(v[1], v[2], v[3]), v[4], v[5])


############################################################
# Movement
############################################################

abstract type Movement end
struct PolarMovement{T} <: Movement
    direction::SVector{3, T}
    dΘ::T
    dΦ::T
end
struct QuatMovement{T} <: Movement
    direction::SVector{3,T}
    rotation::Rotor
end

# convert between movement types
QuatMovement(pm::PolarMovement) = QuatMovement(pm.direction, from_spherical_coordinates(pm.dΘ,pm.dΦ))
function PolarMovement(qm::QuatMovement)
    eulers = to_euler_angles(qv.rotation)
    PolarMovement(qv.direction, -eulers[2], eulers[1] + eulers[3])
end
Base.convert(::Type{QuatMovement}, pm::PolarMovement) = QuatMovement(pm)
Base.convert(::Type{PolarMovement}, qm::QuatMovement) = PolarMovement(qm)

# convert to/from vectors for serialization and interaction with autodiff
to_vector(qv::QuatMovement) = to_vector(PolarMovement(qv))
to_vector(pv::PolarMovement) = [pv.position..., pv.θ, pv.Φ]
from_vector(::Type{QuatMovement}, v) = QuatMovement(from_vector(PolarMovement, v))
from_vector(::Type{PolarMovement}, v) = PolarMovement(SVector(v[1], v[2], v[3]), v[4], v[5])

# attribute access
Base.angle(qm::QuatMovement) = angle(qm.rotation)
Base.angle(pm::PolarMovement) = angle(QuatMovement(pm))
LinearAlgebra.norm(m::Movement) = norm(m.direction)

# apply movement to viewpoint
apply(qv::QuatViewpoint, qm::QuatMovement) = QuatViewpoint(qv.position + qm.direction, qv.orientation * qm.rotation )
apply(qv::QuatViewpoint, pm::PolarMovement) = apply(qv, QuatMovement(pm))

# check if a movement is valid
# TODO: fix rotation without translation
isvalid(m::Movement, d_max=4., rot_max=deg2rad(20))::Bool = norm(m) ≤ d_max && angle(m) ≤ rot_max 

function make_valid(qm::QuatMovement; d_max=4., rot_max=deg2rad(20))
    if norm(qm.direction) > d_max
        direction = qm.direction / norm(qm.direction) * d_max
    else
        direction = qm.direction
    end
    if angle(qm.rotation) > rot_max
        rotation = squad([Rotor(0), qm.rotation], [0., 1.], rot_max / angle(qm))
    else 
        rotation = qm.rotation
    end
    QuatMovement(direction, rotation)
end
make_valid(pm::PolarMovement; d_max=4., rot_max=deg2rad(20)) = PolarMovement(make_valid(QuatMovement(pm), d_max=d_max, rot_max=rot_max))

# interpolation of movement, not sure if useful
function interpolate(qm::QuatMovement, t)::QuatMovement
    rot = squad([Rotor(0), qm.rotation], [0., 1.], t)
    trans = qm.direction * t
    QuatMovement(trans, rot)
end
function interpolate(pm::PolarMovement, t)::PolarMovement
    PolarMovement(interpolate(QuatMovement(pm), t))
end



