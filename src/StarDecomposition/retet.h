#pragma once

#include <tetgen.h>

#include "mesh.h"
#include "volume_mesh.h"

namespace StarDecomposition
{
bool retetrahedrize(VolumeMesh& mesh);
bool retetrahedrize(const Mesh& mesh, VolumeMesh& out);
}