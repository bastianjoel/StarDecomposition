#pragma once

#include <Eigen/Sparse>

#include <random>
#include <queue>

#include "assertion.h"
#include "lp.h"
#include "mesh.h"
#include "sd_boundary_chebyshev.h"
#include "sd_boundary_lp.h"
#include "vectorq.h"
#include "volume_mesh.h"

namespace StarDecomposition
{
std::vector<VolumeMesh> sd(VolumeMesh& N, std::string, int seed = 0);
std::vector<VolumeMesh> sd(Mesh& N, std::string, int seed = 0);

std::vector<Vector3q> decompose(VolumeMesh& mesh, int seed = 0);
}