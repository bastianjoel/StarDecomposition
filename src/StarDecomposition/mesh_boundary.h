#include "OpenMesh/Core/Mesh/Handles.hh"
#include "mesh.h"
#include "vectorq.h"
#include <map>
#include <utility>
#include <vector>

namespace StarDecomposition
{
class MeshBoundary {
public:
    MeshBoundary(Mesh& mesh) : _mesh(mesh) {
        for (auto h : _mesh.boundary_halfedges()) {
            auto ohe = _mesh.opposite_halfedge_handle(h);
            _normals[ohe] = get_normal(ohe);
        }
    }
    ~MeshBoundary() = default;

    int add_boundary_halfedge(OpenMesh::FaceHandle f) {
        int added = 0;
        for (auto he : _mesh.fh_range(f)) {
            auto ohe = _mesh.opposite_halfedge_handle(he);
            if (_normals.find(he) == _normals.end()) {
                _normals[ohe] = get_normal(ohe);
                added++;
            } else {
                _normals.erase(he);
            }
        }

        return added;
    }

    void remove_boundary_halfedge(OpenMesh::FaceHandle f) {
        int removed = 0;
        for (auto he : _mesh.fh_range(f)) {
            auto ohe = _mesh.opposite_halfedge_handle(he);
            if (_normals.find(ohe) == _normals.end()) {
                _normals[he] = get_normal(he);
            } else {
                _normals.erase(ohe);
                removed++;
            }
        }
    }

    const std::vector<OpenMesh::HalfedgeHandle> get_halfedges() {
        std::vector<OpenMesh::HalfedgeHandle> halfedges(_normals.size());
        int i = 0;
        for (auto h : _normals) {
            halfedges[i++] = h.first;
        }
        return halfedges;
    }

    const std::vector<Eigen::Vector3d> get_normals() {
        std::vector<Eigen::Vector3d> normals(_normals.size());
        int i = 0;
        for (auto h : _normals) {
            normals[i++] = _normals[h.first];
        }
        return normals;
    }

    const std::vector<Eigen::Vector3d> get_normals(const Vector3q& fix) {
        std::vector<Eigen::Vector3d> normals(_normals.size());
        int i = 0;
        for (auto h : _normals) {
            normals[i++] = get_normal(h.first, fix);
        }

        return normals;
    }

    const Vector3q& get_fix_vertex() const {
        return _fix_vertex;
    }

    void set_fix_vertex(const Vector3q& fix_vertex) {
        _fix_vertex = fix_vertex;

        for (auto h : _normals) {
            _normals[h.first] = get_normal(h.first);
        }
    }

    const Vector3q get_fix_vertex_normal() {
        Vector3q n = Vector3q::Zero();
        for (auto h : _normals) {
            n += _normals[h.first].cast<mpq_class>() / static_cast<mpq_class>(_normals.size());
        }

        return n;
    }

private:
    Mesh& _mesh;
    Vector3q _fix_vertex;

    std::map<OpenMesh::HalfedgeHandle, Eigen::Vector3d> _normals;

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
}