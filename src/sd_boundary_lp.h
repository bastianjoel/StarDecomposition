#pragma once

#include "lp.h"
#include "mesh.h"
#include "volume_mesh.h"
#include "vectorq.h"
#include <utility>
#ifdef GUI
#include "viewer.h"
#endif


class StarDecompositionBoundaryLp {
public:
  StarDecompositionBoundaryLp(VolumeMesh& m);
  StarDecompositionBoundaryLp(Mesh& m);
  std::vector<Vector3q> centers();
  std::vector<VolumeMesh> components();
#ifdef GUI
  void set_viewer(Viewer* viewer) { _viewer = viewer; }
#endif
private:
#ifdef GUI
  Viewer* _viewer = nullptr;
#endif
  VolumeMesh _originalMesh;

  Mesh _mesh = Mesh();
  bool _computed;
  int _cmpIdx = 0;
  OpenMesh::FPropHandleT<int> _cmp;
  void run();

  std::pair<StarCenterResult, Vector3q> has_valid_center(Mesh& mesh, Vector3q normal);

  Vector3q _cmpNormal = Vector3q();
  Vector3q _cmpCenter = Vector3q();
  Vector3q _cmpFixV = Vector3q();

  Mesh* _currentCmp = nullptr;
  std::vector<Mesh> _cmpMeshes;
  std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _cmpVertexMap;
  std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _meshVertexMap;

  Mesh add_component(const OpenMesh::FaceHandle& startF);
  bool add_face_to_cmp(Mesh& mesh, OpenMesh::FaceHandle& hf);
  void apply_current_component();
};
