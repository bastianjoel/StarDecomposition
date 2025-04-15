#pragma once

#include "mesh.h"
#include "sd_boundary.h"
#include "volume_mesh.h"
#include "vectorq.h"
#include <utility>

class StarDecompositionBoundaryChebyshev : public StarDecompositionBoundary {
public:
    StarDecompositionBoundaryChebyshev(Mesh& m) : StarDecompositionBoundary(m) {}
    StarDecompositionBoundaryChebyshev(VolumeMesh& vMesh) : StarDecompositionBoundary(vMesh) {}

private:
    void init_component(Mesh& mesh, const OpenMesh::FaceHandle& startF);
    int add_face_to_cmp(Mesh& mesh, const OpenMesh::FaceHandle& hf);
    void finalize_component(Mesh& mesh);

private:
    Vector3q _cmpCenter;

    // Global component meshes
    bool move_fix_vertex(Mesh& mesh, MeshBoundary& boundary, bool shrink);
    bool move_vertex_to(Mesh& mesh, MeshBoundary& boundary, const Vector3q& p);
    OpenMesh::FaceHandle get_opposite_face(Mesh& mesh, const OpenMesh::FaceHandle& origin);
    OpenMesh::FaceHandle get_opposite_face(Mesh& mesh, const OpenMesh::VertexHandle& origin);
    OpenMesh::FaceHandle get_opposite_face(Vector3q from, Vector3q normal);
    std::pair<OpenMesh::FaceHandle, Vector3q> get_fix_vertex_pos(Mesh& mesh, const OpenMesh::FaceHandle& hf);
};
