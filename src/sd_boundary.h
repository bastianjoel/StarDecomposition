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
  VolumeMesh& _originalMesh;

  Mesh _mesh = Mesh();
  bool _computed;
  int _cmpIdx = 0;
  Vector3q _cmpCenter;
  OpenMesh::FPropHandleT<Halfface> _originalHf;
  OpenMesh::FPropHandleT<int> _cmp;
  // std::set<Halfface> _surface_hf;
  // std::vector<Halfface> _boundary;
  void run();

  std::pair<bool, Vector3q> has_valid_center(Mesh& mesh);

  // Global component meshes
  std::pair<OpenMesh::FaceHandle, OpenMesh::VertexHandle> _cmpFixV;
  std::vector<Mesh> _cmpMeshes;
  std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _cmpVertexMap;
  Mesh add_component(OpenMesh::FaceHandle startF);
  bool add_face_to_cmp(Mesh& mesh, OpenMesh::FaceHandle& hf);
  bool point_in_triangle(Vector3q p, Vector3q v0, Vector3q v1, Vector3q v2);
  bool triangles_intersect(std::vector<Vector3q> t);
  void apply_current_component();
  OpenMesh::FaceHandle get_opposite_face(Mesh& mesh, OpenMesh::FaceHandle& hf);
  std::pair<OpenMesh::FaceHandle, Vector3q> get_fix_vertex_pos(Mesh& mesh, OpenMesh::FaceHandle& hf);

  Vector3q get_face_center(Mesh& mesh, OpenMesh::FaceHandle& hf);
  // OpenMesh::FaceHandle check_intersecting(OpenMesh::FaceHandle hf);

  // Direct add (experiment)
  bool add_halfface(std::set<Halfface>& cmpHf, Halfface& hf);
};
