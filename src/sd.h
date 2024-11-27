#pragma once

#include <Eigen/Sparse>

#include <random>
#include <queue>

#include "assertion.h"
#include "lp.h"
#include "mesh.h"
#include "sd_boundary.h"
#include "vectorq.h"

std::vector<Mesh> sd(Mesh& N);

std::vector<Vector3q> decompose(Mesh& mesh);
