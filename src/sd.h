#pragma once

#include <Eigen/Sparse>

#include <random>
#include <queue>

#include "assertion.h"
#include "lp.h"
#include "mesh.h"
#include "vectorq.h"

std::vector<Mesh> sdm(Mesh& N);

std::vector<Vector3q> decompose(Mesh& mesh);
