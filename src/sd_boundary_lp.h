#pragma once

#include "lp.h"
#include "mesh.h"
#include "sd_boundary.h"
#include "volume_mesh.h"
#include "vectorq.h"
#include <optional>
#include <utility>

class StarDecompositionBoundaryLp : public StarDecompositionBoundary {
public:
    StarDecompositionBoundaryLp(Mesh& m) : StarDecompositionBoundary(m) {}
    StarDecompositionBoundaryLp(VolumeMesh& vMesh) : StarDecompositionBoundary(vMesh) {}

private:
    std::pair<StarCenterResult, Vector3q> has_valid_center(Mesh& mesh, Vector3q normal);

    Vector3q _cmpNormal = Vector3q();
    Vector3q _cmpCenter = Vector3q();
    Vector3q _cmpFixV = Vector3q();

    Mesh init_component(const OpenMesh::FaceHandle& startF);
    int add_face_to_cmp(Mesh& mesh, const OpenMesh::FaceHandle& hf);
    void finalize_component(Mesh& mesh);

    std::optional<Vector3q> get_fix_vertex_pos(Mesh& mesh, const Vector3q& vPos, const Vector3q& n); 
    bool is_valid_fixpoint(Mesh& mesh, const Vector3q& fixV, Vector3q& center);
    bool is_valid_component(Mesh& mesh, const Vector3q& fixV, Vector3q& center);
};
