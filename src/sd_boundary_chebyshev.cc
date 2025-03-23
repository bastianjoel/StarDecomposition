#include "sd_boundary_chebyshev.h"
#include "assertion.h"
#include <Eigen/src/Core/Matrix.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <OpenMesh/Core/IO/MeshIO.hh>

void StarDecompositionBoundaryChebyshev::finalize_component(Mesh& cmpMesh) {
    Vector3q openingCenter = Vector3q::Zero();
    int numBoundaryVertices = 0;
    for (auto f : cmpMesh.vf_range(_cmpFixV.second)) {
        for (auto h : cmpMesh.fh_range(f)) {
            auto to = cmpMesh.to_vertex_handle(h);
            auto from = cmpMesh.to_vertex_handle(h);
            if (to == _cmpFixV.second || from == _cmpFixV.second) {
                continue;
            }

            openingCenter += cmpMesh.data(to).point_q();
            numBoundaryVertices++;
        }
    }
    openingCenter /= numBoundaryVertices;
    for (int i = 0; i < 3; i++) {
        mpq_class step = i / 4.0;
        auto fixVertexPos = cmpMesh.data(_cmpFixV.second).point_q();
        auto moveTo = openingCenter + step * (fixVertexPos - openingCenter);
        if (move_vertex_to(cmpMesh, _cmpFixV.second, moveTo)) {
            break;
        }
    }

    double maxArea = -1;
    auto fixVertex = _mesh.add_vertex_q(cmpMesh.data(_cmpFixV.second).point_q());
    for (auto f : _mesh.faces()) {
        if (_mesh.property(_selected, f)) {
            _mesh.delete_face(f, true);
            double cMaxArea = _mesh.calc_face_area(f);
            if (cMaxArea > maxArea) {
                maxArea = cMaxArea;
            }
        }
    }

    auto v2 = _mesh.data(fixVertex).point_q();
    for (auto h : _mesh.halfedges()) {
        if (_mesh.is_boundary(h)) {
            auto f = _mesh.add_face(_mesh.from_vertex_handle(h), _mesh.to_vertex_handle(h), fixVertex);
            _mesh.update_normal_q(f);
            _mesh.property(_selected, f) = false;
        }
    }

    _cmpFixV.second = fixVertex;
    cmpMesh.garbage_collection();

    _mesh.delete_isolated_vertices();
    _mesh.garbage_collection();
    _mesh.generate_bvh();
}

bool StarDecompositionBoundaryChebyshev::move_fix_vertex(Mesh& mesh, bool shrink) {
    Vector3q vPos = mesh.data(_cmpFixV.second).point_q();
    auto nD = mesh.calc_normal(_cmpFixV.second) * (shrink ? 1 : -1);
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

        if (move_vertex_to(mesh, _cmpFixV.second, p)) {
            _cmpFixV.first = opposite;
            return true;
        }
    }

    return false;
}

bool StarDecompositionBoundaryChebyshev::move_vertex_to(Mesh& mesh, OpenMesh::VertexHandle& moveV, const Vector3q& p) {
    for (auto vf : mesh.vf_range(moveV)) {
        std::vector<OpenMesh::VertexHandle> vertices;
        for (auto v : mesh.fv_range(vf)) {
            if (moveV != v) {
                vertices.push_back(v);
            }
        }
        std::vector<Vector3q> t = { mesh.data(vertices[0]).point_q(), mesh.data(vertices[1]).point_q(), p };
        if (_mesh.triangle_intersects(t, { _meshVertexMap[vertices[0]], _meshVertexMap[vertices[1]] }).is_valid()) {
            return false;
        }
    }

    Vector3q prevPos = mesh.data(moveV).point_q();
    mesh.move(_cmpFixV.second, p);
    auto cmpCenter = mesh.star_center();
    if (!cmpCenter.first) {
        mesh.move(_cmpFixV.second, prevPos);
        return false;
    }

    _cmpCenter = cmpCenter.second;
    return true;
}

Mesh StarDecompositionBoundaryChebyshev::init_component(const OpenMesh::FaceHandle& startF) {
    Mesh mesh;
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

    auto opposite = get_fix_vertex_pos(_mesh, startF);
    if (opposite.first.is_valid()) {
        auto point = opposite.second;
        auto newVertex = mesh.add_vertex_q(point);
        _cmpFixV = std::make_pair(opposite.first, newVertex);
        for (auto he : mesh.fh_range(face)) {
            auto ohe = mesh.opposite_halfedge_handle(he);
            auto of = mesh.add_face(mesh.from_vertex_handle(ohe), mesh.to_vertex_handle(ohe), newVertex);
            mesh.update_normal_q(of);
        }
    } else {
        std::cout << "Opposite face not found" << std::endl;
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
    return mesh;
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
bool StarDecompositionBoundaryChebyshev::add_face_to_cmp(Mesh& mesh, const OpenMesh::FaceHandle& newFace) {
    // TODO: Check which edge connects to the component
    // TODO: Use connecting edge to remove an temporary face
    // TODO: Connect open boundaries to the component
    // TODO: Check if opposite faces should be expanded
    OpenMesh::VertexHandle newFaceVertex;
    std::vector<OpenMesh::VertexHandle> triangle;
    int deletedHelperFaces = 0;
    OpenMesh::VertexHandle fixedVertex = _cmpFixV.second;
    TxDeleteMesh txMesh(mesh);

    for (auto fh : _mesh.fh_range(newFace)) {
        auto hv = _mesh.from_vertex_handle(fh);
        auto hvTo = _mesh.to_vertex_handle(fh);
        if (_cmpVertexMap[hv].is_valid() && _cmpVertexMap[hvTo].is_valid()) {
            auto h = mesh.find_halfedge(_cmpVertexMap[hv], _cmpVertexMap[hvTo]);
            if (h.is_valid()) {
                auto f = mesh.face_handle(h);
                txMesh.delete_face(f);
                deletedHelperFaces++;
            }
        }
    }

    for (auto fh : _mesh.fh_range(newFace)) {
        auto hv = _mesh.from_vertex_handle(fh);
        if (!_cmpVertexMap[hv].is_valid()) {
            auto newVertex = mesh.add_vertex_q(_mesh.data(hv).point_q());
            _cmpVertexMap[hv] = newVertex;
            _meshVertexMap[newVertex] = hv;
            newFaceVertex = hv;
        }

        triangle.push_back(_cmpVertexMap[hv]);
    }

    bool shouldCheck = false;
    std::vector<OpenMesh::FaceHandle> illegalFaces;
    if (deletedHelperFaces != 0 && (deletedHelperFaces != 1 || newFaceVertex.is_valid())) {
        OpenMesh::FaceHandle face = txMesh.add_face(triangle);
        mesh.set_normal_q(face, _mesh.data(newFace).normal_q());
        for (auto fh : mesh.fh_range(face)) {
            if (mesh.opposite_halfedge_handle(fh).is_boundary()) {
                auto v0 = mesh.to_vertex_handle(fh);
                auto v1 = mesh.from_vertex_handle(fh);
                auto v2 = fixedVertex;
                std::vector<Vector3q> t = { mesh.data(v0).point_q(), mesh.data(v1).point_q(), mesh.data(v2).point_q() };
                if (_mesh.triangle_intersects(t, { _meshVertexMap[v0], _meshVertexMap[v1] }).is_valid()) {
                    illegalFaces.push_back(face);
                }

                auto f = txMesh.add_face(v0, v1, v2);
                mesh.update_normal_q(f);
            }
        }
        shouldCheck = true;
    } else if (deletedHelperFaces == 1 && !newFaceVertex.is_valid()) {
        // std::cout << "Multi add" << std::endl;
        // TODO: Add enclosed faces
        // 1. Determine open edges (2)
        // 2. Check number of enclosed faces + determine enclosed faces
        // 3. Use smaller chunk
        // 4. Check if enclose face can be added => set illegalTriangle
        // 5. Remove enclose faces
        // 6. Add new faces
    }

    bool illegalTriangle = illegalFaces.size() > 0;
    if (illegalTriangle && move_fix_vertex(mesh, true)) {
        illegalTriangle = false;
    }

    std::pair<bool, Vector3q> cmpCenter;
    if (shouldCheck && !illegalTriangle) {
        cmpCenter = mesh.star_center();
        /*
        if (!cmpCenter.first && illegalFaces.size() == 0 && move_fix_vertex(mesh, false)) {
            cmpCenter = std::make_pair(true, _cmpCenter);
        }
        */
    }

#ifdef GUI
    std::vector<Eigen::Vector3d> positions;
    for (auto v : triangle) {
        positions.push_back(mesh.point(v));
    }
#endif

    if (illegalTriangle || !shouldCheck || !cmpCenter.first) {
        if (newFaceVertex.is_valid()) {
            _meshVertexMap.erase(_cmpVertexMap[newFaceVertex]);
            _cmpVertexMap.erase(newFaceVertex);
        }

#ifdef GUI
        if (illegalTriangle) {
            _viewer->add_triangle(positions, { 1, 0, 0 });
        } else if (shouldCheck) {
            _viewer->add_triangle(positions, { 0, 0, 1 });
        } else {
            _viewer->add_triangle(positions, { 1, 0, 1 });
        }
#endif

        txMesh.revert();

        return false;
    } else {
#ifdef GUI
        _viewer->add_triangle(positions, { 0, 1, 0 });
#endif
    }


    _cmpCenter = cmpCenter.second;
    return true;
}

OpenMesh::FaceHandle StarDecompositionBoundaryChebyshev::get_opposite_face(Mesh& mesh, const OpenMesh::FaceHandle& origin) {
    Vector3q vPos = mesh.face_center(origin);
    Vector3q n = -mesh.data(origin).normal_q();

    return get_opposite_face(vPos, n);
}

OpenMesh::FaceHandle StarDecompositionBoundaryChebyshev::get_opposite_face(Mesh& mesh, const OpenMesh::VertexHandle& origin) {
    Vector3q vPos = mesh.data(origin).point_q();
    auto n = mesh.calc_normal(origin);

    return get_opposite_face(vPos, n.cast<mpq_class>());
}

OpenMesh::FaceHandle StarDecompositionBoundaryChebyshev::get_opposite_face(Vector3q vPos, Vector3q n) {
    return _mesh.ray_intersects(vPos, n);
}

std::pair<OpenMesh::FaceHandle, Vector3q> StarDecompositionBoundaryChebyshev::get_fix_vertex_pos(Mesh& mesh, const OpenMesh::FaceHandle& hf) {
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
