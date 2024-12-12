#pragma once

#include "mesh.h"
#include "volume_mesh.h"
#include "vectorq.h"

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
  OpenMesh::FPropHandleT<Halfface> _originalHf;
  OpenMesh::FPropHandleT<int> _cmp;
  OpenMesh::FPropHandleT<Vector3q> _normalQ;
  OpenMesh::VPropHandleT<Vector3q> _Q;
  // std::set<Halfface> _surface_hf;
  // std::vector<Halfface> _boundary;
  void run();

  bool has_valid_center(Mesh mesh);

  // Global component meshes
  std::vector<Mesh> _cmpMeshes;
  std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle> _cmpVertexMap;
  Mesh add_component(OpenMesh::FaceHandle startF);
  bool add_hf_to_cmp(int cmp, OpenMesh::FaceHandle hf);

  // Direct add (experiment)
  bool add_halfface(std::set<Halfface>& cmpHf, Halfface& hf);
};
