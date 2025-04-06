#include "OpenMesh/Core/Mesh/Handles.hh"
#include "mesh.h"
#include "vectorq.h"
#include <map>
#include <set>
#include <vector>

class MeshBoundary {
public:
  MeshBoundary(Mesh& mesh) : _mesh(mesh) {
    for (auto h : _mesh.boundary_halfedges()) {
      auto ohe = _mesh.opposite_halfedge_handle(h);
      _boundaryHalfedges.insert(ohe);
      _normals[ohe] = get_normal(ohe);
    }
  }
  ~MeshBoundary() = default;

  int add_boundary_halfedge(OpenMesh::FaceHandle f) {
    int added = 0;
    for (auto he : _mesh.fh_range(f)) {
      auto ohe = _mesh.opposite_halfedge_handle(he);
      if (_boundaryHalfedges.find(he) == _boundaryHalfedges.end()) {
        _boundaryHalfedges.insert(ohe);
        _normals[ohe] = get_normal(ohe);
        added++;
      } else {
        _boundaryHalfedges.erase(he);
        _normals.erase(he);
      }
    }

    return added;
  }

  void remove_boundary_halfedge(OpenMesh::FaceHandle f) {
    int removed = 0;
    for (auto he : _mesh.fh_range(f)) {
      auto ohe = _mesh.opposite_halfedge_handle(he);
      if (_boundaryHalfedges.find(ohe) == _boundaryHalfedges.end()) {
        _boundaryHalfedges.insert(he);
        _normals[he] = get_normal(he);
      } else {
        _boundaryHalfedges.erase(ohe);
        _normals.erase(ohe);
        removed++;
      }
    }
  }

  const std::vector<Eigen::Vector3d> get_normals() {
    std::vector<Eigen::Vector3d> normals(_boundaryHalfedges.size());
    int i = 0;
    for (auto h : _boundaryHalfedges) {
      normals[i++] = _normals[h];
    }
    return normals;
  }

  const std::vector<Eigen::Vector3d> get_normals(const Vector3q& fix) {
    std::vector<Eigen::Vector3d> normals(_boundaryHalfedges.size());
    int i = 0;
    for (auto h : _boundaryHalfedges) {
      normals[i++] = get_normal(h, fix);
    }

    return normals;
  }

  const Vector3q& get_fix_vertex() const {
    return _fix_vertex;
  }

  void set_fix_vertex(const Vector3q& fix_vertex) {
    _fix_vertex = fix_vertex;

    for (auto h : _boundaryHalfedges) {
        _normals[h] = get_normal(h);
    }
  }

private:
  Mesh& _mesh;
  Vector3q _fix_vertex;

  std::map<OpenMesh::HalfedgeHandle, Eigen::Vector3d> _normals;
  std::set<OpenMesh::HalfedgeHandle> _boundaryHalfedges; 

  Eigen::Vector3d get_normal(const OpenMesh::HalfedgeHandle& h) {
    return get_normal(h, _fix_vertex);
  }

  Eigen::Vector3d get_normal(const OpenMesh::HalfedgeHandle& h, const Vector3q& fix_vertex) {
    Vector3q v0 = _mesh.data(_mesh.to_vertex_handle(h)).point_q();
    Vector3q v1 = _mesh.data(_mesh.from_vertex_handle(h)).point_q();
    Eigen::Vector3d n = (v1 - v0).cross(fix_vertex - v0).unaryExpr([](mpq_class x) { return x.get_d(); });
    return (-n).normalized();
  }
};
