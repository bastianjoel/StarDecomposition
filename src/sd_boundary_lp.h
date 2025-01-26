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
  bool _wasVolumeMesh;
  VolumeMesh _originalMesh;

  Mesh _mesh = Mesh();
  bool _computed;
  int _cmpIdx = 0;
  OpenMesh::FPropHandleT<Halfface> _originalHf;
  OpenMesh::FPropHandleT<int> _cmp;
  // std::set<Halfface> _surface_hf;
  // std::vector<Halfface> _boundary;
  void run();

  std::pair<bool, Vector3q> has_valid_center(Mesh& mesh);
  std::pair<StarCenterResult, Vector3q> has_valid_center(Mesh& mesh, const Vector3q& point, std::vector<std::pair<Vector3q, Vector3q>> additional = std::vector<std::pair<Vector3q, Vector3q>>());

  std::pair<OpenMesh::VertexHandle, Vector3q> _cmpFixV = std::make_pair(OpenMesh::VertexHandle(-1), Vector3q());
  std::vector<Mesh> _cmpMeshes;
  Mesh add_component(const OpenMesh::FaceHandle& startF);

  // Global component meshes
  Mesh* _currentCmp = nullptr;
  std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _cmpVertexMap;
  std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _meshVertexMap;
  bool add_face_to_cmp(Mesh& mesh, OpenMesh::FaceHandle& hf);
  void apply_current_component();

  OpenMesh::FaceHandle triangles_intersect(std::vector<Vector3q> t, std::vector<OpenMesh::VertexHandle> v);

  // OpenMesh::FaceHandle check_intersecting(OpenMesh::FaceHandle hf);

  // Direct add (experiment)
  bool add_halfface(std::set<Halfface>& cmpHf, Halfface& hf);
};
