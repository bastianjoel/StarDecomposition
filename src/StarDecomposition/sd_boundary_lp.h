#pragma once

#include "lp.h"
#include "mesh.h"
#include "sd_boundary.h"
#include "vectorq.h"
#include "volume_mesh.h"
#include <optional>
#include <utility>

namespace StarDecomposition
{
class StarDecompositionBoundaryLp : public StarDecompositionBoundary {
public:
    StarDecompositionBoundaryLp(Mesh& m, int seed = 0) : StarDecompositionBoundary(m, seed) {}
    StarDecompositionBoundaryLp(VolumeMesh& vMesh, int seed = 0) : StarDecompositionBoundary(vMesh, seed) {}

private:
    std::pair<StarCenterResult, Vector3q> has_valid_center(Mesh& mesh, Vector3q normal);

    Vector3q _cmpNormal = Vector3q();
    Vector3q _cmpCenter = Vector3q();

    void init_component(Mesh& mesh, const OpenMesh::FaceHandle& startF);
    int add_face_to_cmp(Mesh& mesh, const OpenMesh::FaceHandle& hf);
    void finalize_component(Mesh& mesh);

    std::optional<Vector3q> get_fix_vertex_pos(Mesh& mesh, MeshBoundary &boundary, const Vector3q& vPos, const Vector3q& n); 
    bool is_valid_fixpoint(Mesh& mesh, MeshBoundary &boundary, const Vector3q& fixV, Vector3q& center);
    bool is_valid_component(Mesh& mesh, MeshBoundary &boundary, const Vector3q& fixV, Vector3q& center);
};
}