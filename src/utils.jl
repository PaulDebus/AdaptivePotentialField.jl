# SPDX-FileCopyrightText: 2023 Paul Debus <paul.debus@uni-weimar.de>
#
# SPDX-License-Identifier: GPL-3.0-or-later

export random_points_and_normals, loadmesh

using Meshes

function random_points_and_normals(rng::AbstractRNG, Ω::Domain,
    method::HomogeneousSampling)
	size    = method.size
	weights = measure.(Ω)

	# sample elements with weights proportial to measure
	w = WeightedSampling(size, weights, replace=true)

	# within each element sample a single point
	h = HomogeneousSampling(1)

	tris = sample(rng, Ω, w)

	(first(sample(rng, e, h)) for e in tris) , tris
end

function loadmesh(path)
    convert(Meshes.SimpleMesh, FileIO.load(path))
end

function Base.convert(::Type{Meshes.SimpleMesh}, mesh::GeometryBasics.Mesh)
    # load a mesh as a Meshes.jl mesh, normal loading loads as GeometryBasics mesh
	points = [Tuple(p) for p in Set(mesh.position)]
	indices = Dict(p => i for (i, p) in enumerate(points))
	connectivities = map(mesh) do el
		Meshes.connect(Tuple(indices[Tuple(p)] for p in el))
	end
	Meshes.SimpleMesh(Meshes.Point.(points), connectivities)
end
