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
    bool _computed = false;
    std::vector<VolumeMesh> _components;
    std::vector<Vector3q> _centers;

    Mesh _mesh = Mesh();
    OpenMesh::FPropHandleT<int> _cmp;

    virtual void run() = 0;
    virtual Mesh add_component(const OpenMesh::FaceHandle& startF) = 0;
    virtual bool add_face_to_cmp(Mesh& mesh, const OpenMesh::FaceHandle& hf) = 0;

#ifdef GUI
protected:
    Viewer* _viewer = nullptr;
public:
    void set_viewer(Viewer* viewer) { _viewer = viewer; }
#endif
};
