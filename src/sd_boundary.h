#pragma once

#include "mesh.h"
#include "volume_mesh.h"
#include "vectorq.h"
#include <utility>

class StarDecompositionBoundary {
public:
  StarDecompositionBoundary(VolumeMesh& m);
  StarDecompositionBoundary(Mesh& m);
  std::vector<Vector3q> centers();
  std::vector<VolumeMesh> components();
private:
  bool _wasVolumeMesh;
  VolumeMesh _originalMesh;

  int noCenter = 0;
  int triangleIntersect = 0;
  int other = 0;
  Mesh _mesh = Mesh();
  bool _computed;
  int _cmpIdx = 0;
  Vector3q _cmpCenter;
  OpenMesh::FPropHandleT<Halfface> _originalHf;
  OpenMesh::FPropHandleT<int> _cmp;
  // std::set<Halfface> _surface_hf;
  // std::vector<Halfface> _boundary;
  void run();

  // Global component meshes
  std::pair<OpenMesh::FaceHandle, OpenMesh::VertexHandle> _cmpFixV;
  std::vector<Mesh> _cmpMeshes;
  Mesh* _currentCmp = nullptr;
  std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _cmpVertexMap;
  std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _meshVertexMap;
  Mesh add_component(OpenMesh::FaceHandle startF);
  bool add_face_to_cmp(Mesh& mesh, OpenMesh::FaceHandle& hf);
  bool move_fix_vertex(Mesh& mesh);
  bool move_vertex_to(Mesh& mesh, OpenMesh::VertexHandle& v, const Vector3q& p);
  void apply_current_component();
  OpenMesh::FaceHandle get_opposite_face(Mesh& mesh, OpenMesh::FaceHandle& origin);
  OpenMesh::FaceHandle get_opposite_face(Mesh& mesh, OpenMesh::VertexHandle& origin);
  OpenMesh::FaceHandle get_opposite_face(Vector3q from, Vector3q normal);
  std::pair<OpenMesh::FaceHandle, Vector3q> get_fix_vertex_pos(Mesh& mesh, OpenMesh::FaceHandle& hf);

  // OpenMesh::FaceHandle check_intersecting(OpenMesh::FaceHandle hf);

  // Direct add (experiment)
  bool add_halfface(std::set<Halfface>& cmpHf, Halfface& hf);
};
