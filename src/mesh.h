#pragma once

#include <Eigen/Dense>
#include "OpenMesh/Core/Mesh/Attributes.hh"
#include "OpenMesh/Core/Mesh/Traits.hh"
#include <Eigen/src/Core/Matrix.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Geometry/EigenVectorT.hh>
#include "assertion.h"
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
public:
    Vector3q face_center(OpenMesh::FaceHandle fh) {
        Vector3q center = Vector3q(0, 0, 0);
        for (auto v : fv_range(fh)) {
            center += data(v).point_q();
        }
        return center / 3;
    }

    // https://math.stackexchange.com/questions/51326/determining-if-an-arbitrary-point-lies-inside-a-triangle-defined-by-three-points
    // https://gdbooks.gitbooks.io/3dcollisions/content/Chapter4/point_in_triangle.html
    bool point_on_face(OpenMesh::FaceHandle fh, Vector3q p) {
        std::vector<Vector3q> t;
        for (auto hv : fv_range(fh)) {
            t.push_back(data(hv).point_q());
        }
        t[0] -= p;
        t[1] -= p;
        t[2] -= p;

        auto u = t[1].cross(t[2]);
        auto v = t[2].cross(t[0]);
        auto w = t[0].cross(t[1]);

        if (u.dot(v) < 0 || u.dot(w) < 0) {
            return false;
        }

        return true;
    }

    OpenMesh::SmartVertexHandle add_vertex_q(const Vector3q p) {
        auto fixVertex = add_vertex(p.unaryExpr([](mpq_class x) { return x.get_d(); }));
        data(fixVertex).set_point_q(p);

        return fixVertex;
    }

    Vector3q update_normal_q(OpenMesh::FaceHandle fh) {
        std::vector<Vector3q> triangle;
        for (auto fv : fv_range(fh)) {
            triangle.push_back(data(fv).point_q());
        }
        Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
        set_normal_q(fh, n);

        return n;
    }

    void set_normal_q(OpenMesh::FaceHandle fh, const Vector3q& n) {
        data(fh).set_normal_q(n);
        set_normal(fh, n.unaryExpr([](mpq_class x) { return x.get_d(); }));
    }

    void move(OpenMesh::VertexHandle vh, const Vector3q& p) {
        data(vh).set_point_q(p);
        set_point(vh, p.unaryExpr([](mpq_class x) { return x.get_d(); }));
        for (auto f : vf_range(vh)) {
            update_normal_q(f);
        }
    }

    OpenMesh::FaceHandle triangle_intersects(std::vector<Vector3q> t, std::vector<OpenMesh::VertexHandle> borderVertices) {
        Vector3q n = (t[1] - t[0]).cross(t[2] - t[0]);
        for (auto face : this->faces()) {
            if (triangle_intersects(t, n, borderVertices, face)) {
                return face;
            }
        }

        return OpenMesh::FaceHandle();
    }

    bool triangle_intersects(std::vector<Vector3q> t, Vector3q n, std::vector<OpenMesh::VertexHandle> borderVertices, OpenMesh::FaceHandle face) {
        std::vector<OpenMesh::VertexHandle> fv = {
            this->to_vertex_handle(this->halfedge_handle(face)),
            this->to_vertex_handle(this->next_halfedge_handle(this->halfedge_handle(face))),
            this->to_vertex_handle(this->next_halfedge_handle(this->next_halfedge_handle(this->halfedge_handle(face))))
        };

        // Move connecting vertices slightly to avoid intersection
        std::vector<bool> shareV = { false, false, false };
        for (auto v : borderVertices) {
            for (int i = 0; i < 3; i++) {
                if (v == fv[i]) {
                    shareV[i] = true;
                }
            }
        }
        short sharedVertices = shareV[0] + shareV[1] + shareV[2];
        ASSERT(sharedVertices < 3);

        Vector3q fNormal = this->data(face).normal_q();
        std::vector<Vector3q> vq = { this->data(fv[0]).point_q(), this->data(fv[1]).point_q(), this->data(fv[2]).point_q() };
        if (shareV[0] && shareV[1]) {
            vq[0] = vq[0] + (vq[2] - vq[0]) * 1e-6;
            vq[1] = vq[1] + (vq[2] - vq[1]) * 1e-6;
        } else if (shareV[0] && shareV[2]) {
            vq[0] = vq[0] + (vq[1] - vq[0]) * 1e-6;
            vq[2] = vq[2] + (vq[1] - vq[2]) * 1e-6;
        } else if (shareV[1] && shareV[2]) {
            vq[1] = vq[1] + (vq[0] - vq[1]) * 1e-6;
            vq[2] = vq[2] + (vq[0] - vq[2]) * 1e-6;
        } else if (shareV[0]) {
            vq[0] = vq[0] + (vq[2] - vq[0]) * 1e-6 + (vq[1] - vq[0]) * 1e-6;
        } else if (shareV[1]) {
            vq[1] = vq[1] + (vq[2] - vq[1]) * 1e-6 + (vq[0] - vq[1]) * 1e-6;
        } else if (shareV[2]) {
            vq[2] = vq[2] + (vq[1] - vq[2]) * 1e-6 + (vq[0] - vq[2]) * 1e-6;
        }

        return tri_tri_intersect(t[0], t[1], t[2], n, vq[0], vq[1], vq[2], fNormal);
    }

    std::pair<bool, Vector3q> star_center() {
        std::vector<Eigen::Vector3d> positions;
        std::vector<Eigen::Vector3d> normals;
        for (auto face : this->faces()) {
            if (face.is_valid() && !face.deleted()) {
                positions.push_back(this->point(this->to_vertex_handle(this->halfedge_handle(face))));
                normals.push_back(this->normal(face).normalized());
            }
        }

        Vector3q newCenter = kernel_chebyshev_center(positions, normals).cast<mpq_class>();
        for (auto face : this->faces()) {
            std::vector<Vector3q> triangle;
            for (auto fv : this->fv_range(face)) {
                triangle.push_back(this->data(fv).point_q());
            }
            Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
            if ((triangle[0] - newCenter).dot(n) <= 0) {
                return std::make_pair(false, Vector3q());
            }
        }

        return std::make_pair(true, newCenter);
    }

    std::pair<bool, mpq_class> intersection_factor(const Vector3q& p, const Vector3q& dir, OpenMesh::FaceHandle face) {
        auto faceNormal = this->data(face).normal_q();
        auto div = faceNormal.transpose() * dir;
        if (div[0] == 0) {
            return std::make_pair(false, mpq_class());
        }

        // Calculate intersection
        auto c = faceNormal.transpose() * this->data(this->to_vertex_handle(this->halfedge_handle(face))).point_q();
        auto nominator = faceNormal.transpose() * p - c;
        auto denominator = faceNormal.transpose() * dir;
        return std::make_pair(true, (-nominator / denominator)[0]);
    }
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
