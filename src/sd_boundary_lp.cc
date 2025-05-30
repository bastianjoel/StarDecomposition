#include "sd_boundary_lp.h"
#include "OpenMesh/Core/Mesh/Handles.hh"
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

    for (auto b : _boundaries) {
        auto cmpFixVertex = cmpMesh.add_vertex_q(b.get_fix_vertex());
        for (auto h : b.get_halfedges()) {
            cmpMesh.add_face(cmpMesh.from_vertex_handle(h), cmpMesh.to_vertex_handle(h), cmpFixVertex);
        }

        auto fixVertex = _mesh.add_vertex_q(b.get_fix_vertex());
        for (auto h : b.get_halfedges()) {
            auto fromV = _meshVertexMap[cmpMesh.to_vertex_handle(h)];
            auto toV = _meshVertexMap[cmpMesh.from_vertex_handle(h)];
            auto f = _mesh.add_face(fromV, toV, fixVertex);
            _mesh.update_normal_q(f);
            _mesh.property(_selected, f) = false;
            _mesh.property(_origBound, f) = false;
        }
    }

    cmpMesh.garbage_collection();

    _mesh.delete_isolated_vertices();
    _mesh.garbage_collection();
    _mesh.generate_bvh();
}

void StarDecompositionBoundaryLp::init_component(Mesh& mesh, const OpenMesh::FaceHandle& startF) {
    std::vector<OpenMesh::VertexHandle> newHfVertices;
    for (auto hv : _mesh.fv_range(startF)) {
        auto newVertex = mesh.add_vertex_q(_mesh.data(hv).point_q());
        _cmpVertexMap[hv] = newVertex;
        _meshVertexMap[newVertex] = hv;
        newHfVertices.push_back(newVertex);
    }

    auto face = mesh.add_face(newHfVertices);
    mesh.set_normal_q(face, _mesh.data(startF).normal_q());
    _boundaries[0].add_boundary_halfedge(face);

    _cmpNormal = _mesh.data(startF).normal_q();
    _cmpCenter = _mesh.face_center(startF);
    auto p = get_fix_vertex_pos(mesh, _boundaries[0], _cmpCenter, _cmpNormal);
    if (p.has_value()) {
        _boundaries[0].set_fix_vertex(p.value());
    } else {
        _boundaries[0].set_fix_vertex(_cmpCenter - _mesh.data(startF).normal_q() * 1e-4);
    }

    _mesh.property(_selected, startF) = true;
}

/**
 * Adds a face to an existing triangle mesh
 */
int StarDecompositionBoundaryLp::add_face_to_cmp(Mesh& mesh, const OpenMesh::FaceHandle& newFace) {
    OpenMesh::VertexHandle newFaceVertex;
    std::vector<OpenMesh::VertexHandle> triangle;

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

    MeshBoundary& boundary = _boundaries[0];
    bool valid = is_valid_with(boundary, newFace, _cmpCenter);
    Eigen::Vector3d color = { (double)!valid, (double)valid, 0 };
    auto addedFacesNormal = _mesh.data(newFace).normal_q();

    TxDeleteMesh txMesh(mesh);
    OpenMesh::FaceHandle face = txMesh.add_face(triangle);
    mesh.set_normal_q(face, _mesh.data(newFace).normal_q());
    boundary.add_boundary_halfedge(face);
    mpq_class f1 = _nextComponentFaces / (_nextComponentFaces + 1.0f);
    mpq_class f2 = 1.0f / (_nextComponentFaces + 1.0f);
    Vector3q meshNormal = ((_cmpNormal * f1) + (addedFacesNormal * f2));// .unaryExpr([](mpq_class x) { return x.get_d(); }).cast<mpq_class>();

    bool invalidCenter = false;
    if (!valid) {
        auto meshNormalD = meshNormal.unaryExpr([](mpq_class x) { return x.get_d(); }).cast<mpq_class>();
        auto cmpCenter = has_valid_center(mesh, meshNormalD);
        if (cmpCenter.first == INVALID) {
            invalidCenter = true;
            color = { 0, 0, 1 };
        } else {
            auto p = get_fix_vertex_pos(mesh, boundary, cmpCenter.second, meshNormalD);
            if (p.has_value()) {
                valid = true;
                // TODO: Verbosely set in is_valid_component
                // _cmpCenter = cmpCenter.second;
                boundary.set_fix_vertex(p.value());
                _recheckFailed = true;
                color = { 0, 1, 1 };
            }
        }
    }

    if (!valid) {
        if (newFaceVertex.is_valid()) {
            _meshVertexMap.erase(_cmpVertexMap[newFaceVertex]);
            _cmpVertexMap.erase(newFaceVertex);
        }

        boundary.remove_boundary_halfedge(face);
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

std::optional<Vector3q> StarDecompositionBoundaryLp::get_fix_vertex_pos(Mesh& mesh, MeshBoundary &boundary, const Vector3q& cPos, const Vector3q& n) {
    mpq_class t;
    auto normal = n;
    auto center = cPos;
    auto opposite = _mesh.ray_intersects(cPos, -normal, t);
    if (!opposite.is_valid()) {
        return std::nullopt;
    }

    if ((_mesh.data(*_mesh.fv_begin(opposite)).point_q() - _cmpCenter).dot(_mesh.data(opposite).normal_q()) > 0) {
        Vector3q p;
        for (int i = 0; i < 4; i++) {
            bool valid = true;
            p = cPos - normal * (t / pow(2, i));
            p = p.unaryExpr([](mpq_class x) { return x.get_d(); }).cast<mpq_class>();

            if (is_valid_fixpoint(mesh, boundary, p, center)) {
                _cmpCenter = center.cast<mpq_class>();
                return p;
            }
        }
    }

    return std::nullopt;
}

bool StarDecompositionBoundaryLp::is_valid_fixpoint(Mesh& mesh, MeshBoundary &boundary, const Vector3q& fixV, Vector3q& center) {
    if (!is_valid_component(mesh, boundary, fixV, center)) {
        return false;
    }

    for (auto b : _boundaries) {
        for (auto h : b.get_halfedges()) {
            if (triangle_intersects(h, fixV)) {
                return false;
            }
        }
    }

    return true;
}

std::pair<StarCenterResult, Vector3q> StarDecompositionBoundaryLp::has_valid_center(Mesh& mesh, Vector3q normal) {
    std::vector<Eigen::Vector3d> positions(_nextComponentFaces + 1);
    std::vector<Eigen::Vector3d> normals(_nextComponentFaces + 1);
    int i = 0;
    for (auto face : mesh.faces()) {
        positions[i] = mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(face)));
        normals[i] = mesh.data(face).normal_normalized();
        i++;
    }

    auto result = star_center_close_to(normal.unaryExpr([](mpq_class x) { return x.get_d(); }).normalized(), positions, normals);
    if (result.first == INVALID) {
        return std::make_pair(result.first, Vector3q::Zero());
    }

    Vector3q newCenter = result.second.cast<mpq_class>();
    return std::make_pair(result.first, newCenter);
}

bool StarDecompositionBoundaryLp::is_valid_component(Mesh& mesh, MeshBoundary &boundary, const Vector3q& fixV, Vector3q& center) {
    std::vector<Eigen::Vector3d> positions(_nextComponentFaces + 1);
    std::vector<Eigen::Vector3d> normals(_nextComponentFaces + 1);
    int i = 0;
    for (auto face : mesh.faces()) {
        positions[i] = mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(face)));
        normals[i] = mesh.data(face).normal_normalized();
        i++;
    }

    for (auto b : _boundaries) {
        for (auto n : b.get_normals(fixV)) {
            positions.push_back(fixV.unaryExpr([](mpq_class x) { return x.get_d(); }));
            normals.push_back(n);
        }
    }

    Eigen::Vector3d centerD = center.unaryExpr([](mpq_class x) { return x.get_d(); });
    bool valid = true;
    for (int i = positions.size() - 1; i >= 0; i--) {
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

    center = newCenter.cast<mpq_class>();
    return true;
}
