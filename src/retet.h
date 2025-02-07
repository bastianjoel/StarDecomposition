#pragma once

#include <tetgen.h>

#include "mesh.h"
#include "volume_mesh.h"

bool retetrahedrize(VolumeMesh& mesh);
bool retetrahedrize(Mesh& mesh, VolumeMesh& out);
