#pragma once

#include "mesh.h"
#include "mesh_boundary.h"
#include "volume_mesh.h"
#include "vectorq.h"

#ifdef GUI
#include "viewer.h"
#endif

class StarDecompositionBoundary {
public:
    StarDecompositionBoundary(VolumeMesh& m);
    StarDecompositionBoundary(Mesh& m);

    std::vector<Vector3q> centers();
    std::vector<VolumeMesh> components();

protected:
    // When set to true recoverable faces are rechecked starting from next iteration
    bool _recheckFailed = false;

    // Result of the decomposition
    std::vector<VolumeMesh> _components;
    std::vector<Vector3q> _centers;

    // Mesh to decompose
    Mesh _mesh = Mesh();
    OpenMesh::FPropHandleT<int> _selected;

    // Current processing component
    Mesh _nextComponent = Mesh();
    int _nextComponentFaces = 0;
    std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _cmpVertexMap;
    std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _meshVertexMap;
    std::vector<MeshBoundary> _boundaries;

    // Add components to the result
    void add_component(const Mesh& mesh);
    void fallback(const Mesh& mesh);

    // Algorithm specific methods
    virtual void init_component(Mesh& mesh, const OpenMesh::FaceHandle& startF) = 0;
    virtual int add_face_to_cmp(Mesh& mesh, const OpenMesh::FaceHandle& hf) = 0;
    virtual void finalize_component(Mesh& mesh) = 0;

private:
    bool _computed = false;
    void start();


#ifdef GUI
protected:
    Viewer* _viewer = nullptr;
public:
    void set_viewer(Viewer* viewer) { _viewer = viewer; }
#endif
};
