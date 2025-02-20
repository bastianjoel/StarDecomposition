#pragma once

#include "mesh.h"
#include "viewer.h"
#include "volume_mesh.h"
#include "vectorq.h"

class StarDecompositionBoundary {
public:
    StarDecompositionBoundary(VolumeMesh& m);
    StarDecompositionBoundary(Mesh& m);

    std::vector<Vector3q> centers();
    std::vector<VolumeMesh> components();

protected:
    // Result of the decomposition
    std::vector<VolumeMesh> _components;
    std::vector<Vector3q> _centers;

    // Mesh to decompose
    Mesh _mesh = Mesh();
    OpenMesh::FPropHandleT<int> _selected;

    // Current processing component
    Mesh _nextComponent = Mesh();
    std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _cmpVertexMap;
    std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _meshVertexMap;

    // Add components to the result
    void add_component(const Mesh& mesh);
    void fallback(const Mesh& mesh);

    // Algorithm specific methods
    virtual Mesh init_component(const OpenMesh::FaceHandle& startF) = 0;
    virtual bool add_face_to_cmp(Mesh& mesh, const OpenMesh::FaceHandle& hf) = 0;
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
