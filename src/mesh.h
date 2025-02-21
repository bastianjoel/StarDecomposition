#pragma once

#include <Eigen/Dense>
#include "OpenMesh/Core/Mesh/Attributes.hh"
#include "OpenMesh/Core/Mesh/Traits.hh"
#include <Eigen/src/Core/Matrix.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/EigenVectorT.hh>
#include "assertion.h"
#include "bvh.h"
#include "lp.h"
#include "vectorq.h"
#include "tritri.h"

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

// typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  Mesh;
class Mesh : public OpenMesh::TriMesh_ArrayKernelT<MyTraits> {
private:
    OpenMesh::HalfedgeHandle lastBoundaryHalfedge;
    BVHNode* _bvh = nullptr;
    OpenMesh::FaceHandle triangle_intersects_bvh(BVHNode* node, const std::vector<Vector3q>& t, const Vector3q& n, const std::vector<OpenMesh::VertexHandle>& borderVertices);
    OpenMesh::FaceHandle ray_intersects_bvh(BVHNode* node, const Vector3q& p, const Vector3q& n);
public:
    ~Mesh() {
        if (_bvh != nullptr) {
            // delete _bvh;
        }
    }

    Vector3q face_center(OpenMesh::FaceHandle fh);

    OpenMesh::SmartVertexHandle add_vertex_q(const Vector3q p);
    Vector3q update_normal_q(OpenMesh::FaceHandle fh);
    void set_normal_q(OpenMesh::FaceHandle fh, const Vector3q& n);

    void move(OpenMesh::VertexHandle vh, const Vector3q& p);
    std::pair<bool, Vector3q> star_center();
    std::pair<bool, mpq_class> intersection_factor(const Vector3q& p, const Vector3q& dir, OpenMesh::FaceHandle face);
    std::vector<OpenMesh::HalfedgeHandle> boundary_halfedges();

    bool point_on_face(OpenMesh::FaceHandle fh, Vector3q p);
    OpenMesh::FaceHandle get_face_in_dir(const Vector3q& vPos, const Vector3q& n);

    OpenMesh::FaceHandle triangle_intersects(const std::vector<Vector3q>& t, const std::vector<OpenMesh::VertexHandle>& borderVertices);
    bool triangle_intersects(const std::vector<Vector3q>& t, const Vector3q& n, const std::vector<OpenMesh::VertexHandle>& borderVertices, const OpenMesh::FaceHandle& face);

    OpenMesh::FaceHandle ray_intersects(const Vector3q& o, const Vector3q& n);

    void generate_bvh();
    BVHNode* generate_bvh(std::vector<OpenMesh::FaceHandle>& faces, int depth = 0);
};

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
