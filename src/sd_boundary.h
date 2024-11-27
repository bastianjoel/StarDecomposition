#pragma once

#include "mesh.h"
#include "vectorq.h"

class StarDecompositionBoundary {
public:
  StarDecompositionBoundary(Mesh& m);
  std::vector<Vector3q> centers();
  std::vector<Mesh> components();
private:
  Mesh& _mesh;
  bool _computed;
  Property<Halfface, int> _cmp;
  Property<Halfface, Vector3q> _normal;
  std::set<Halfface> _surface_hf;
  void run();
  bool can_add_halfface(Mesh& cmpMesh, Halfface& hf, Property<Vertex, Vertex> vertexMap);
};
