#pragma once

#include "mesh.h"
#include "volume_mesh.h"
#include "vectorq.h"

class StarDecompositionBoundary {
public:
  StarDecompositionBoundary(VolumeMesh& m);
  std::vector<Vector3q> centers();
  std::vector<VolumeMesh> components();
private:
  VolumeMesh& _mesh;
  bool _computed;
  int _cmpIdx = 0;
  Property<Vertex, Vector3q> _Q;
  Property<Halfface, Vector3q> _normal;
  Property<Halfface, int> _cmp;
  std::set<Halfface> _surface_hf;
  std::vector<Halfface> _boundary;
  void run();

  bool has_valid_center(MyMesh mesh);

  // Global component meshes
  std::vector<MyMesh> _cmpMeshes;
  std::map<Vertex, OpenMesh::SmartVertexHandle> _cmpVertexMap;
  MyMesh add_component(Halfface& startHf);
  bool add_hf_to_cmp(int cmp, Halfface& hf);

  // Direct add (experiment)
  bool add_halfface(std::set<Halfface>& cmpHf, Halfface& hf);
};
