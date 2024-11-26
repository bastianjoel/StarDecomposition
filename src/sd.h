#pragma once

#include <Eigen/Sparse>

#include <random>
#include <queue>

#include "assertion.h"
#include "lp.h"
#include "mesh.h"
#include "vectorq.h"

std::vector<Mesh> sd(Mesh& N);

std::vector<Vector3q> decompose(Mesh& mesh);

std::vector<Vector3q> decompose_by_corners(Mesh& mesh);

std::set<Halfface> get_surface_halffaces(Mesh& mesh);

bool can_add_halfface(Mesh& parentMesh, Mesh& cmpMesh, Halfface& hf, Property<Vertex, Vertex> vertexMap);
