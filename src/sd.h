#pragma once

#include <Eigen/Sparse>

#include <random>
#include <queue>

#include "assertion.h"
#include "lp.h"
#include "volume_mesh.h"
#include "sd_boundary.h"
#include "vectorq.h"

std::vector<VolumeMesh> sd(VolumeMesh& N);

std::vector<Vector3q> decompose(VolumeMesh& mesh);
