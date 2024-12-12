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

  VertexAttributes(OpenMesh::Attributes::Status);

  VertexTraits {
    private:
      Vector3q _q;
    public:
      const Vector3q& Q() const { return _q; }
      void set_q(const Vector3q& q) { _q = q; }
  };
};

typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  Mesh;
