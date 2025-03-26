#include "sd_boundary_lp.h"
#include "OpenMesh/Core/Mesh/Handles.hh"
#include "assertion.h"
#include "lp.h"
#include "vectorq.h"
#include <Eigen/src/Core/Matrix.h>
#include <cmath>
#include <cstdio>
#include <map>
#include <OpenMesh/Core/IO/MeshIO.hh>

void StarDecompositionBoundaryLp::finalize_component(Mesh& cmpMesh) {
    for (auto f : _mesh.faces()) {
        if (_mesh.property(_selected, f)) {
            _mesh.delete_face(f, true);
        }
    }

    auto cmpFixVertex = cmpMesh.add_vertex_q(_cmpFixV);
    for (auto h : cmpMesh.boundary_halfedges()) {
        cmpMesh.add_face(cmpMesh.from_vertex_handle(h), cmpMesh.to_vertex_handle(h), cmpFixVertex);
    }

    auto fixVertex = _mesh.add_vertex_q(_cmpFixV);
    for (auto h : _mesh.boundary_halfedges()) {
        auto f = _mesh.add_face(_mesh.from_vertex_handle(h), _mesh.to_vertex_handle(h), fixVertex);
        _mesh.update_normal_q(f);
        _mesh.property(_selected, f) = false;
    }

    cmpMesh.garbage_collection();

    _mesh.delete_isolated_vertices();
    _mesh.garbage_collection();
    _mesh.generate_bvh();
}

/*
std::pair<OpenMesh::FaceHandle, Vector3q> StarDecompositionBoundaryLp::get_fix_vertex_pos(Mesh& mesh, const OpenMesh::FaceHandle& hf) {
    Vector3q vPos = mesh.face_center(hf);
    Vector3q n = -mesh.data(hf).normal_q();

    auto opposite = _mesh.ray_intersects(vPos, n);
    if (!opposite.is_valid()) {
        return std::make_pair(opposite, Vector3q::Zero());
    }

    n = mesh.data(hf).normal_q();
    mpq_class r = mesh.intersection_factor(vPos, n, opposite).second / 2;
    Vector3q p = vPos + r * n;

    bool intersects;
    do {
        intersects = false;
        for (auto fh : mesh.fh_range(hf)) {
            auto v0 = mesh.to_vertex_handle(fh);
            auto v1 = mesh.from_vertex_handle(fh);
            std::vector<Vector3q> t = { mesh.data(v0).point_q(), mesh.data(v1).point_q(), p };
            if (_mesh.triangle_intersects(t, { v0, v1 }).is_valid()) {
                intersects = true;
                r /= 2;
                p = vPos + r * n;
                ASSERT(r > 1e-6 || r < -1e-6);
                break;
            }
        }
    } while (intersects);

    // TODO: Check if there is a way to make this workaround unnecessary
    p = p.unaryExpr([](mpq_class x) { return x.get_d(); }).cast<mpq_class>();

    return std::make_pair(opposite, p);
}
*/

Mesh StarDecompositionBoundaryLp::init_component(const OpenMesh::FaceHandle& startF) {
    Mesh mesh;
    std::vector<OpenMesh::VertexHandle> newHfVertices;
    for (auto hv : _mesh.fv_range(startF)) {
        auto newVertex = mesh.add_vertex_q(_mesh.data(hv).point_q());
        _cmpVertexMap[hv] = newVertex;
        _meshVertexMap[newVertex] = hv;
        newHfVertices.push_back(newVertex);
    }

    auto face = mesh.add_face(newHfVertices);
    mesh.set_normal_q(face, _mesh.data(startF).normal_q());

    _cmpNormal = _mesh.data(startF).normal_q();
    _cmpCenter = _mesh.face_center(startF);
    auto p = get_fix_vertex_pos(mesh, _cmpCenter, _cmpNormal);
    if (p.has_value()) {
        _cmpFixV = p.value();
    } else {
        _cmpFixV = _cmpCenter - _mesh.data(startF).normal_q() * 1e-4;
    }

    _mesh.property(_selected, startF) = true;
    return mesh;
}

/**
 * Adds a face to an existing triangle mesh
 */
int StarDecompositionBoundaryLp::add_face_to_cmp(Mesh& mesh, const OpenMesh::FaceHandle& newFace) {
    OpenMesh::VertexHandle newFaceVertex;
    std::vector<OpenMesh::VertexHandle> triangle;
    TxDeleteMesh txMesh(mesh);

    OpenMesh::HalfedgeHandle singleExistingEdge;
    for (auto fh : _mesh.fh_range(newFace)) {
        auto hv = _mesh.from_vertex_handle(fh);
        auto hvTo = _mesh.to_vertex_handle(fh);
        if (!_cmpVertexMap[hv].is_valid()) {
            auto newVertex = mesh.add_vertex_q(_mesh.data(hv).point_q());
            _cmpVertexMap[hv] = newVertex;
            _meshVertexMap[newVertex] = hv;
            newFaceVertex = hv;
        } else {
            OpenMesh::HalfedgeHandle h = mesh.find_halfedge(_cmpVertexMap[hvTo], _cmpVertexMap[hv]);
            if (h.is_valid()) {
                if (singleExistingEdge.is_valid()) {
                    singleExistingEdge = OpenMesh::HalfedgeHandle();
                } else {
                    singleExistingEdge = h;
                }
            }
        }

        triangle.push_back(_cmpVertexMap[hv]);
    }

    // Face would split mesh into two components
    if (singleExistingEdge.is_valid() && !newFaceVertex.is_valid()) {
        return 1;
    }

    Eigen::Vector3d color = { 1, 0, 0 };
    bool invalidCenter = false;
    mpq_class addedFaces = 1;
    auto addedFacesNormal = _mesh.data(newFace).normal_q();
    OpenMesh::FaceHandle face = txMesh.add_face(triangle);
    mesh.set_normal_q(face, _mesh.data(newFace).normal_q());
    addedFacesNormal /= addedFaces;

    mpq_class f1 = (mesh.n_faces() - addedFaces) / mesh.n_faces();
    mpq_class f2 = addedFaces / mesh.n_faces();
    Vector3q meshNormal = ((_cmpNormal * f1) + (addedFacesNormal * f2));// .unaryExpr([](mpq_class x) { return x.get_d(); }).cast<mpq_class>();

    bool valid = true;
    if ((_mesh.data(*_mesh.fv_begin(newFace)).point_q() - _cmpCenter).dot(_mesh.data(newFace).normal_q()) > 0) {
        for (auto fh : mesh.fh_range(face)) {
            if (mesh.opposite_halfedge_handle(fh).is_boundary()) {
                auto v0 = mesh.to_vertex_handle(fh);
                auto v1 = mesh.from_vertex_handle(fh);
                auto n = (mesh.data(v1).point_q() - mesh.data(v0).point_q()).cross(_cmpFixV - mesh.data(v0).point_q());
                if ((mesh.data(v0).point_q() - _cmpCenter).dot(n) <= 0) {
                    color = { 0.5, 0, 0 };
                    valid = false;
                    break;
                }

                std::vector<Vector3q> t = { mesh.data(v0).point_q(), mesh.data(v1).point_q(), _cmpFixV };
                if (_mesh.triangle_intersects(t, { _meshVertexMap[v0], _meshVertexMap[v1] }).is_valid()) {
                    color = { 1, 0, 0 };
                    valid = false;
                    break;
                }
            }
        }

        if (valid) {
            color = { 0, 1, 0 };
        }
    } else {
        valid = false;
    }

    if (!valid) {
        auto cmpCenter = has_valid_center(mesh, meshNormal);
        if (cmpCenter.first == INVALID) {
            invalidCenter = true;
            color = { 0, 0, 1 };
        } else {
            auto p = get_fix_vertex_pos(mesh, cmpCenter.second, meshNormal);
            if (p.has_value()) {
                valid = true;
                // TODO: Verbosely set in is_valid_component
                // _cmpCenter = cmpCenter.second;
                _cmpFixV = p.value();
                _recheckFailed = true;
                color = { 0, 1, 1 };
#ifdef GUI
                // _viewer->clear_extras();
                _viewer->add_sphere(_cmpCenter.unaryExpr([](mpq_class x) { return x.get_d(); }), 0.001, { 0, 1, 0 });
                _viewer->add_sphere(_cmpFixV.unaryExpr([](mpq_class x) { return x.get_d(); }), 0.001, { 0, 1, 1 });
#endif
            }
        }
    }

    if (!valid) {
        if (newFaceVertex.is_valid()) {
            _meshVertexMap.erase(_cmpVertexMap[newFaceVertex]);
            _cmpVertexMap.erase(newFaceVertex);
        }

        txMesh.revert();
    } else {
        _cmpNormal = meshNormal;
    }

#ifdef GUI
    std::vector<Eigen::Vector3d> positions;
    for (auto v : triangle) {
        positions.push_back(mesh.point(v));
    }
    _viewer->add_triangle(positions, color);
#endif

    if (valid) {
        return 0;
    } else if (invalidCenter) {
        return 3;
    } else {
        return 2;
    }
}

std::optional<Vector3q> StarDecompositionBoundaryLp::get_fix_vertex_pos(Mesh& mesh, const Vector3q& cPos, const Vector3q& n) {
    mpq_class t;
    auto normal = n;
    auto opposite = _mesh.ray_intersects(cPos, -normal, t);
    if ((_mesh.data(*_mesh.fv_begin(opposite)).point_q() - _cmpCenter).dot(_mesh.data(opposite).normal_q()) > 0) {
        Vector3q p;
        auto boundary = mesh.boundary_halfedges();
        for (int i = 0; i < 4; i++) {
            bool valid = true;
            p = cPos - normal * (t / pow(2, i));
            p = p.unaryExpr([](mpq_class x) { return x.get_d(); }).cast<mpq_class>();

            for (auto h : boundary) {
                auto v0 = mesh.to_vertex_handle(h);
                auto v1 = mesh.from_vertex_handle(h);
                std::vector<Vector3q> t = { mesh.data(v0).point_q(), mesh.data(v1).point_q(), p };
                auto intersects = _mesh.triangle_intersects(t, { _meshVertexMap[v0], _meshVertexMap[v1] });
                if (intersects.is_valid()) {
                    valid = false;
                    break;
                }
            }

            if (valid && is_valid_component(mesh, p, cPos)) {
                return p;
            }
        }
    }

    return std::nullopt;
}

std::pair<StarCenterResult, Vector3q> StarDecompositionBoundaryLp::has_valid_center(Mesh& mesh, Vector3q normal) {
    std::vector<Eigen::Vector3d> positions;
    std::vector<Eigen::Vector3d> normals;
    for (auto face : mesh.faces()) {
        if (face.is_valid() && !face.deleted()) {
            positions.push_back(mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(face))));
            normals.push_back(mesh.normal(face).normalized());
        }
    }

    auto result = star_center_close_to(normal.unaryExpr([](mpq_class x) { return x.get_d(); }).normalized(), positions, normals);
    if (result.first == INVALID) {
        return std::make_pair(result.first, Vector3q::Zero());
    }

    Vector3q newCenter = result.second.cast<mpq_class>();
    return std::make_pair(result.first, newCenter);
}

bool StarDecompositionBoundaryLp::is_valid_component(Mesh& mesh, const Vector3q& fixV, const Vector3q& center) {
    std::vector<Eigen::Vector3d> positions;
    std::vector<Eigen::Vector3d> normals;
    for (auto face : mesh.faces()) {
        if (face.is_valid() && !face.deleted()) {
            positions.push_back(mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(face))));
            normals.push_back(mesh.normal(face).normalized());
        }
    }

    for (auto h : mesh.boundary_halfedges()) {
        positions.push_back(fixV.unaryExpr([](mpq_class x) { return x.get_d(); }));
        Vector3q v0 = mesh.data(mesh.to_vertex_handle(h)).point_q();
        Vector3q v1 = mesh.data(mesh.from_vertex_handle(h)).point_q();
        Eigen::Vector3d n = (v1 - v0).cross(fixV - v0).unaryExpr([](mpq_class x) { return x.get_d(); });
        normals.push_back((-n).normalized());
    }

    Eigen::Vector3d centerD = center.unaryExpr([](mpq_class x) { return x.get_d(); });
    bool valid = true;
    for (int i = 0; i < positions.size(); i++) {
        if ((positions[i] - centerD).dot(normals[i]) <= 0) {
            valid = false;
            break;
        }
    }

    if (valid) {
        _cmpCenter = center;
        return true;
    }

    auto newCenter = kernel_chebyshev_center(positions, normals);
    for (int i = 0; i < positions.size(); i++) {
        if ((positions[i] - newCenter).dot(normals[i]) <= 0) {
            return false;
        }
    }

    _cmpCenter = newCenter.cast<mpq_class>();
    return true;
}
