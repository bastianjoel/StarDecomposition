#include "sd_boundary_chebyshev.h"
#include "assertion.h"
#include <Eigen/src/Core/Matrix.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <optional>

void StarDecomposition::StarDecompositionBoundaryChebyshev::finalize_component(Mesh& cmpMesh) {
    for (auto b : _boundaries) {
        Vector3q openingCenter = Vector3q::Zero();
        int numBoundaryVertices = 0;
        for (auto h : b.get_halfedges()) {
            auto to = cmpMesh.to_vertex_handle(h);
            openingCenter += cmpMesh.data(to).point_q();
            numBoundaryVertices++;
        }
        openingCenter /= numBoundaryVertices;
        for (int i = 0; i < 3; i++) {
            mpq_class step = i / 4.0;
            auto fixVertexPos = b.get_fix_vertex();
            auto moveTo = openingCenter + step * (fixVertexPos - openingCenter);
            if (move_vertex_to(cmpMesh, b, moveTo)) {
                break;
            }
        }
    }

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

void StarDecomposition::StarDecompositionBoundaryChebyshev::init_component(Mesh& mesh, const OpenMesh::FaceHandle& startF) {
    std::vector<OpenMesh::VertexHandle> newHfVertices;
    for (auto hv : _mesh.fv_range(startF)) {
        auto hvQ = _mesh.data(hv).point_q();
        auto newVertex = mesh.add_vertex_q(hvQ);
        _cmpVertexMap[hv] = newVertex;
        _meshVertexMap[newVertex] = hv;
        newHfVertices.push_back(newVertex);
    }

    auto face = mesh.add_face(newHfVertices);
    mesh.set_normal_q(face, _mesh.data(startF).normal_q());
    _boundaries[0].add_boundary_halfedge(face);

    auto opposite = get_fix_vertex_pos(_mesh, startF);
    if (opposite.first.is_valid()) {
        _boundaries[0].set_fix_vertex(opposite.second);
    } else {
        std::cerr << "Opposite face not found" << std::endl;
        {
            if (!OpenMesh::IO::write_mesh(_mesh, "debug/_error.obj")) {
                std::cerr << "write error\n";
                exit(1);
            }
        }
        exit(1);
    }

#ifdef GUI
    std::vector<Eigen::Vector3d> positions;
    for (auto v : _mesh.fv_range(startF)) {
        positions.push_back(_mesh.point(v));
    }
    _viewer->add_triangle(positions, { 0, 1, 0 });
#endif

    _mesh.property(_selected, startF) = true;
}

/**
 * Adds a face to an existing triangle mesh
 *
 * Handles the following cases:
 * 1. Face introduces a new vertex
 *  In this case two new faces will be added to the mesh
 * 2. Face connects to an existing vertex
 *  In this case the face will be added to the mesh
 */
int StarDecomposition::StarDecompositionBoundaryChebyshev::add_face_to_cmp(Mesh& mesh, const OpenMesh::FaceHandle& newFace) {
    // TODO: Check which edge connects to the component
    // TODO: Use connecting edge to remove an temporary face
    // TODO: Connect open boundaries to the component
    // TODO: Check if opposite faces should be expanded
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
        } else if (_cmpVertexMap[hvTo].is_valid()) {
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
    bool illegalTriangle = false;
    if (!valid) {
        for (auto fh : _mesh.fh_range(newFace)) {
            auto hv = _mesh.from_vertex_handle(fh);
            auto hvTo = _mesh.to_vertex_handle(fh);
            if (!_nextComponent.find_halfedge(_cmpVertexMap[hv], _cmpVertexMap[hvTo]).is_valid()) {
                auto v0 = _mesh.data(hvTo).point_q();
                auto v1 = _mesh.data(hv).point_q();
                std::vector<Vector3q> t = { v0, v1, boundary.get_fix_vertex() };
                if (_mesh.triangle_intersects(t, { hv, hvTo }).is_valid()) {
                    illegalTriangle = true;
                    break;
                }
            }
        }
    }

    TxDeleteMesh txMesh(mesh);
    OpenMesh::FaceHandle face = txMesh.add_face(triangle);
    mesh.set_normal_q(face, _mesh.data(newFace).normal_q());
    std::vector<OpenMesh::FaceHandle> illegalFaces;
    boundary.add_boundary_halfedge(face);
    Vector3q originalFixVertex = boundary.get_fix_vertex();

    if (illegalTriangle && move_fix_vertex(mesh, boundary, true)) {
        illegalTriangle = false;
    }

    std::optional<Vector3q> cmpCenter = std::nullopt;
    if (valid) {
        cmpCenter = _cmpCenter;
    } else if (!illegalTriangle) {
        cmpCenter = is_next_component_valid();
    }

#ifdef GUI
    std::vector<Eigen::Vector3d> positions;
    for (auto v : triangle) {
        positions.push_back(mesh.point(v));
    }
#endif

    if (!cmpCenter.has_value()) {
        if (newFaceVertex.is_valid()) {
            _meshVertexMap.erase(_cmpVertexMap[newFaceVertex]);
            _cmpVertexMap.erase(newFaceVertex);
        }

#ifdef GUI
        if (illegalTriangle) {
            _viewer->add_triangle(positions, { 1, 0, 0 });
        } else {
            _viewer->add_triangle(positions, { 1, 0, 1 });
        }
#endif

        boundary.remove_boundary_halfedge(face);
        boundary.set_fix_vertex(originalFixVertex);
        txMesh.revert();

        return 1;
    } else {
#ifdef GUI
        _viewer->add_triangle(positions, { 0, 1, 0 });
#endif
    }

    _cmpCenter = cmpCenter.value();
    return 0;
}

bool StarDecomposition::StarDecompositionBoundaryChebyshev::move_fix_vertex(Mesh& mesh, MeshBoundary& boundary, bool shrink) {
    Vector3q vPos = boundary.get_fix_vertex();
    auto nD = boundary.get_fix_vertex_normal() * (shrink ? 1 : -1); // TODO: Check if direction correct
    Vector3q n = nD.cast<mpq_class>();

    auto opposite = get_opposite_face(vPos, -n);
    if (!opposite.is_valid()) {
        std::cout << "Opposite face not found" << std::endl;
        return false;
    }

    auto r = _mesh.intersection_factor(vPos, n, opposite);
    if (!r.first) {
        return false;
    }

    r.second /= 4;

    for (int i = 1; i < 4; i++) {
        Vector3q p =  vPos + r.second * i * n;
        p = p.unaryExpr([](mpq_class x) { return x.get_d(); }).cast<mpq_class>();

        if (move_vertex_to(mesh, boundary, p)) {
            return true;
        }
    }

    return false;
}

bool StarDecomposition::StarDecompositionBoundaryChebyshev::move_vertex_to(Mesh& mesh, MeshBoundary& boundary, const Vector3q& p) {
    for (auto h : boundary.get_halfedges()) {
        if (triangle_intersects(h, p)) {
            return false;
        }
    }

    Vector3q prevPos = boundary.get_fix_vertex();
    boundary.set_fix_vertex(p);
    auto cmpCenter = is_next_component_valid();
    if (!cmpCenter.has_value()) {
        boundary.set_fix_vertex(prevPos);
        return false;
    }

    _cmpCenter = cmpCenter.value();
    return true;
}

OpenMesh::FaceHandle StarDecomposition::StarDecompositionBoundaryChebyshev::get_opposite_face(Mesh& mesh, const OpenMesh::FaceHandle& origin) {
    Vector3q vPos = mesh.face_center(origin);
    Vector3q n = -mesh.data(origin).normal_q();

    return get_opposite_face(vPos, n);
}

OpenMesh::FaceHandle StarDecomposition::StarDecompositionBoundaryChebyshev::get_opposite_face(Mesh& mesh, const OpenMesh::VertexHandle& origin) {
    Vector3q vPos = mesh.data(origin).point_q();
    auto n = mesh.calc_normal(origin);

    return get_opposite_face(vPos, n.cast<mpq_class>());
}

OpenMesh::FaceHandle StarDecomposition::StarDecompositionBoundaryChebyshev::get_opposite_face(Vector3q vPos, Vector3q n) {
    return _mesh.ray_intersects(vPos, n);
}

std::pair<OpenMesh::FaceHandle, Vector3q> StarDecomposition::StarDecompositionBoundaryChebyshev::get_fix_vertex_pos(Mesh& mesh, const OpenMesh::FaceHandle& hf) {
    auto opposite = get_opposite_face(mesh, hf);
    if (!opposite.is_valid()) {
        return std::make_pair(opposite, Vector3q::Zero());
    }

    Vector3q vPos = mesh.face_center(hf);
    auto n = mesh.data(hf).normal_q();
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
