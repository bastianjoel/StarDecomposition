#pragma once

#include "lp.h"
#include "mesh.h"
#include "volume_mesh.h"
#include "vectorq.h"
#include <utility>

class StarDecompositionBoundaryLp {
public:
  StarDecompositionBoundaryLp(VolumeMesh& m);
  StarDecompositionBoundaryLp(Mesh& m);
  std::vector<Vector3q> centers();
  std::vector<VolumeMesh> components();
private:
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
  std::vector<Mesh> _cmpMeshes;
  Mesh add_component(const OpenMesh::FaceHandle& startF);

  // Global component meshes
  Mesh* _currentCmp = nullptr;
  std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _cmpVertexMap;
  std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _meshVertexMap;
  bool add_face_to_cmp(Mesh& mesh, OpenMesh::FaceHandle& hf);
  void apply_current_component();
};
