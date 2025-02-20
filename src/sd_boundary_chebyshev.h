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
    void run();
    Mesh add_component(const OpenMesh::FaceHandle& startF);
    bool add_face_to_cmp(Mesh& mesh, const OpenMesh::FaceHandle& hf);

private:
    int noCenter = 0;
    int triangleIntersect = 0;
    int other = 0;
    int _cmpIdx = 0;
    Vector3q _cmpCenter;

    // Global component meshes
    std::pair<OpenMesh::FaceHandle, OpenMesh::VertexHandle> _cmpFixV;
    std::vector<Mesh> _cmpMeshes;
    Mesh* _currentCmp = nullptr;
    std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _cmpVertexMap;
    std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _meshVertexMap;
    bool move_fix_vertex(Mesh& mesh, bool shrink);
    bool move_vertex_to(Mesh& mesh, OpenMesh::VertexHandle& v, const Vector3q& p);
    void apply_current_component();
    OpenMesh::FaceHandle get_opposite_face(Mesh& mesh, const OpenMesh::FaceHandle& origin);
    OpenMesh::FaceHandle get_opposite_face(Mesh& mesh, const OpenMesh::VertexHandle& origin);
    OpenMesh::FaceHandle get_opposite_face(Vector3q from, Vector3q normal);
    std::pair<OpenMesh::FaceHandle, Vector3q> get_fix_vertex_pos(Mesh& mesh, const OpenMesh::FaceHandle& hf);
};
