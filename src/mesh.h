#pragma once

#include "OpenMesh/Core/Mesh/Attributes.hh"
#include "OpenMesh/Core/Mesh/Traits.hh"
#include "OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh"
#include "vectorq.h"
#include <Eigen/src/Core/Matrix.h>

// define traits
struct MyTraits : public OpenMesh::DefaultTraits
{
  typedef Eigen::Vector3d Point;
  typedef Eigen::Vector3d Normal;

  // use face normals
  FaceAttributes(OpenMesh::Attributes::Normal | OpenMesh::Attributes::Status);

  HalfedgeAttributes(OpenMesh::Attributes::Status);

  EdgeAttributes(OpenMesh::Attributes::Status);

  VertexAttributes(OpenMesh::Attributes::Status);

  VertexTraits {
    private:
      Vector3q _point_q;
    public:
      const Vector3q& point_q() const { return _point_q; }
      void set_point_q(const Vector3q& q) { _point_q = q; }
  };

  FaceTraits {
    private:
      Vector3q _normal_q;
    public:
      const Vector3q& normal_q() const { return _normal_q; }
      void set_normal_q(const Vector3q& q) { _normal_q = q; }
  };
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  Mesh;

template <class Mesh> class TxDeleteT
{
private:
  Mesh& _mesh;
  std::vector<OpenMesh::FaceHandle> _addedFaces = std::vector<OpenMesh::FaceHandle>();
  std::vector<OpenMesh::VertexHandle> _deletedFaces = std::vector<OpenMesh::VertexHandle>();
public:
  // construct with a given mesh
  TxDeleteT(Mesh& mesh) : _mesh(mesh) {}

  OpenMesh::FaceHandle add_face(std::vector<OpenMesh::VertexHandle>& vhandles) {
    return add_face(vhandles[0], vhandles[1], vhandles[2]);
  }

  OpenMesh::FaceHandle add_face(OpenMesh::VertexHandle v0, OpenMesh::VertexHandle v1, OpenMesh::VertexHandle v2) {
    auto face = _mesh.add_face(v0, v1, v2);
    _addedFaces.push_back(face);
    return face;
  }

  void delete_face(OpenMesh::FaceHandle fh) {
    std::vector<OpenMesh::VertexHandle> dTri;
    for (auto v : _mesh.fv_range(fh)) {
      _deletedFaces.push_back(v);
    }
    _mesh.delete_face(fh, true);
  }

  void revert() {
    for (auto f : _addedFaces) {
        _mesh.delete_face(f, true);
    }

    for (int i = 0; i < _deletedFaces.size(); i += 3) {
        auto nF = _mesh.add_face(_deletedFaces[i], _deletedFaces[i + 1], _deletedFaces[i + 2]);
        _mesh.set_normal(nF, _mesh.calc_normal(nF));
    }
  }
};

typedef TxDeleteT<Mesh>  TxDeleteMesh;
