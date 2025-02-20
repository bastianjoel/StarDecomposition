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
  bool add_face_to_cmp(Mesh& mesh, const OpenMesh::FaceHandle& hf);
  void apply_current_component();
};
