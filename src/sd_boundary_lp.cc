#include "sd_boundary_lp.h"
#include "OpenMesh/Core/Mesh/Handles.hh"
#include "lp.h"
#include "retet.h"
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
    for (auto h : cmpMesh.halfedges()) {
        if (cmpMesh.is_boundary(h)) {
            cmpMesh.add_face(cmpMesh.from_vertex_handle(h), cmpMesh.to_vertex_handle(h), cmpFixVertex);
        }
    }

    auto fixVertex = _mesh.add_vertex_q(_cmpFixV);
    for (auto h : _mesh.halfedges()) {
        if (_mesh.is_boundary(h)) {
            auto f = _mesh.add_face(_mesh.from_vertex_handle(h), _mesh.to_vertex_handle(h), fixVertex);
            _mesh.update_normal_q(f);
            _mesh.property(_selected, f) = false;
        }
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
    _cmpFixV = _cmpCenter - _mesh.data(startF).normal_q() * 1e-4;

    _mesh.property(_selected, startF) = true;
    return mesh;
}

/**
 * Adds a face to an existing triangle mesh
 */
bool StarDecompositionBoundaryLp::add_face_to_cmp(Mesh& mesh, const OpenMesh::FaceHandle& newFace) {
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
    std::vector<OpenMesh::FaceHandle> lineFaces;
    if (singleExistingEdge.is_valid() && !newFaceVertex.is_valid()) {
        // Check if one direction consists of a "line" of faces
        auto vExFrom = _meshVertexMap[mesh.from_vertex_handle(singleExistingEdge)];
        auto vExTo = _meshVertexMap[mesh.to_vertex_handle(singleExistingEdge)];
        auto hStart = _mesh.find_halfedge(vExTo, vExFrom);
        auto fStart = newFace;
        bool isLine = true;
        for (int i = 0; i < 2; i++) {
            lineFaces.clear();
            isLine = true;
            hStart = _mesh.next_halfedge_handle(hStart);
            OpenMesh::FaceHandle prev = _mesh.face_handle(hStart);
            OpenMesh::FaceHandle curr = _mesh.opposite_face_handle(hStart);
            OpenMesh::FaceHandle next;
            do {
                int openNeighbours = 0;
                next = OpenMesh::FaceHandle();
                for (auto fh : _mesh.fh_range(curr)) {
                    if (!_mesh.property(_selected, _mesh.opposite_face_handle(fh))) {
                        openNeighbours++;

                        if (_mesh.opposite_face_handle(fh) != prev) {
                            next = _mesh.opposite_face_handle(fh);
                        }
                    }
                }

                if (openNeighbours > 2) {
                    isLine = false;
                    break;
                }

                lineFaces.push_back(curr);
                prev = curr;
                curr = next;
            } while (next.is_valid());

            if (isLine) {
                break;
            }
        }

        if (!isLine) {
            return false;
        }
        // return false;
    }

    bool valid = true;
    bool invalidCenter = false;
    int addedFaces = 1;
    OpenMesh::FaceHandle face = txMesh.add_face(triangle);
    mesh.set_normal_q(face, _mesh.data(newFace).normal_q());
    for (auto f : lineFaces) {
        std::vector<OpenMesh::VertexHandle> lineTriangle;
        for (auto fh : _mesh.fh_range(f)) {
            auto hv = _mesh.from_vertex_handle(fh);
            lineTriangle.push_back(_cmpVertexMap[hv]);
        }

        auto newFace = txMesh.add_face(lineTriangle);
        mesh.set_normal_q(newFace, _mesh.data(f).normal_q());
        addedFaces++;
    }

    Vector3q meshNormal = (_cmpNormal * (mesh.n_faces() - addedFaces) / mesh.n_faces()) + (_mesh.data(newFace).normal_q() * (addedFaces / mesh.n_faces()));
    if ((_mesh.data(*_mesh.fv_begin(newFace)).point_q() - _cmpCenter).dot(_mesh.data(newFace).normal_q()) <= 0) {
LABEL:
        auto cmpCenter = has_valid_center(mesh, meshNormal);
        if (cmpCenter.first == INVALID) {
            invalidCenter = true;
            valid = false;
        } else {
            // TODO: Center needs to be moved
            auto boundary = mesh.boundary_halfedges();
            mpq_class t;
            auto opposite = _mesh.ray_intersects(cmpCenter.second, -meshNormal, t);

            Vector3q p;
            for (int i = 1; i < 4; i++) {
                valid = true;
                p = cmpCenter.second - meshNormal * (t / pow(2, i));
                for (auto h : boundary) {
                    auto v0 = mesh.to_vertex_handle(h);
                    auto v1 = mesh.from_vertex_handle(h);
                    std::vector<Vector3q> t = { mesh.data(v0).point_q(), mesh.data(v1).point_q(), p };
                    auto intersects = _mesh.triangle_intersects(t, { _meshVertexMap[v0], _meshVertexMap[v1] });
                    if (intersects.is_valid()) {
                        cmpCenter.second = p;
                        valid = false;
                        break;
                    }
                }

                if (valid) {
                    break;
                }
            }

            if (valid) {
                _cmpCenter = cmpCenter.second;
                _cmpFixV = p;
#ifdef GUI
                // _viewer->clear_extras();
                _viewer->add_sphere(_cmpCenter.unaryExpr([](mpq_class x) { return x.get_d(); }), 0.01, { 0, 1, 0 });
                _viewer->add_sphere(_cmpFixV.unaryExpr([](mpq_class x) { return x.get_d(); }), 0.01, { 0, 1, 1 });
#endif
            }
        }
    } else {
        for (auto fh : mesh.fh_range(face)) {
            if (mesh.opposite_halfedge_handle(fh).is_boundary()) {
                auto v0 = mesh.to_vertex_handle(fh);
                auto v1 = mesh.from_vertex_handle(fh);
                std::vector<Vector3q> t = { mesh.data(v0).point_q(), mesh.data(v1).point_q(), _cmpFixV };
                if (_mesh.triangle_intersects(t, { _meshVertexMap[v0], _meshVertexMap[v1] }).is_valid()) {
                    goto LABEL;
                    valid = false;
                    break;
                }
            }
        }
    }

#ifdef GUI
    std::vector<Eigen::Vector3d> positions;
    for (auto v : triangle) {
        positions.push_back(mesh.point(v));
    }
#endif

    if (valid == false) {
        if (newFaceVertex.is_valid()) {
            _meshVertexMap.erase(_cmpVertexMap[newFaceVertex]);
            _cmpVertexMap.erase(newFaceVertex);
        }

        txMesh.revert();
#ifdef GUI
        if (invalidCenter) {
            _viewer->add_triangle(positions, { 0, 0, 1 });
        } else {
            _viewer->add_triangle(positions, { 1, 0, 0 });
        }
#endif
    } else {
#ifdef GUI
        _viewer->add_triangle(positions, { 0, 1, 0 });
        for (auto f : lineFaces) {
            std::vector<Eigen::Vector3d> lPos;
            for (auto fh : _mesh.fh_range(f)) {
                lPos.push_back(_mesh.point(_mesh.from_vertex_handle(fh)));
            }
            _viewer->add_triangle(lPos, { 0, 1, 0 });
        }
#endif
        for (auto f : lineFaces) {
            _mesh.property(_selected, f) = true;
        }
    }

    _cmpNormal = meshNormal;
    return valid;
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
