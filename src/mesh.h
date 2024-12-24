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
