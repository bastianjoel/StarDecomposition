#pragma once

#include "lp.h"
#include "mesh.h"
#include "sd_boundary.h"
#include "volume_mesh.h"
#include "vectorq.h"
#include <utility>

class StarDecompositionBoundaryLp : public StarDecompositionBoundary {
public:
    StarDecompositionBoundaryLp(Mesh& m) : StarDecompositionBoundary(m) {}
    StarDecompositionBoundaryLp(VolumeMesh& vMesh) : StarDecompositionBoundary(vMesh) {}

private:
    std::pair<StarCenterResult, Vector3q> has_valid_center(Mesh& mesh, Vector3q& normal);

    Vector3q _cmpNormal = Vector3q();
    Vector3q _cmpCenter = Vector3q();
    Vector3q _cmpFixV = Vector3q();

    Mesh init_component(const OpenMesh::FaceHandle& startF);
    bool add_face_to_cmp(Mesh& mesh, const OpenMesh::FaceHandle& hf);
    void finalize_component(Mesh& mesh);
};
