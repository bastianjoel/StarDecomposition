#pragma once

#include <Eigen/Sparse>

#include <random>
#include <queue>

#include "assertion.h"
#include "lp.h"
#include "volume_mesh.h"
#include "sd_boundary_chebyshev.h"
#include "sd_boundary_lp.h"
#include "vectorq.h"
#include "mesh.h"

std::vector<VolumeMesh> sd(VolumeMesh& N, std::string);
std::vector<VolumeMesh> sd(Mesh& N, std::string);

std::vector<Vector3q> decompose(VolumeMesh& mesh);
