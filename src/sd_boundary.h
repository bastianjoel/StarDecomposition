#pragma once

#include "mesh.h"
#include "vectorq.h"

// define traits
struct MyTraits : public OpenMesh::DefaultTraits
{
  typedef Eigen::Vector3d Point;
  typedef Eigen::Vector3d Normal;

  // use face normals
  FaceAttributes(OpenMesh::Attributes::Normal);

  VertexTraits {
    private:
      Vector3q _q;
    public:
      const Vector3q& Q() const { return _q; }
      void set_q(const Vector3q& q) { _q = q; }
  };
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  MyMesh;

class StarDecompositionBoundary {
public:
  StarDecompositionBoundary(Mesh& m);
  std::vector<Vector3q> centers();
  std::vector<Mesh> components();
private:
  Mesh& _mesh;
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
