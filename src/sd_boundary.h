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
  Property<Vertex, Vector3q> _Q;
  Property<Halfface, Vector3q> _normal;
  Property<Halfface, int> _cmp;
  std::set<Halfface> _surface_hf;
  std::vector<Halfface> _boundary;
  void run();
  bool add_halfface(std::set<Halfface>& cmpHf, Halfface& hf);
};
