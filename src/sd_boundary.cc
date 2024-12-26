#include "sd_boundary.h"
#include "lp.h"
#include "tritri.h"
#include <Eigen/src/Core/Matrix.h>
#include <cmath>
#include <cstdio>
#include <map>
#include <queue>
#include <stack>
#include <OpenMesh/Core/IO/MeshIO.hh>

const int cmpNotSetIdx = -1;
StarDecompositionBoundary::StarDecompositionBoundary(VolumeMesh& mesh) : _computed(false), _wasVolumeMesh(true), _originalMesh(mesh) {
    _mesh.add_property(_cmp);

    _mesh.add_property(_originalHf);

    auto Q = mesh.property<Vertex, Vector3q>("Q");
    std::map<Vertex, OpenMesh::VertexHandle> vMap;
    for (auto v : mesh.vertices()) {
        vMap[v] = _mesh.add_vertex(Mesh::Point(Q[v][0].get_d(), Q[v][1].get_d(), Q[v][2].get_d()));
        _mesh.data(vMap[v]).set_point_q(Q[v]);
    }

    for (auto hf : mesh.boundary_halffaces()) {
        std::vector<OpenMesh::VertexHandle> newHfVertices;
        std::vector<Vector3q> triangle;
        for (auto hv : mesh.halfface_vertices(hf)) {
            newHfVertices.push_back(vMap[hv]);
            triangle.push_back(Q[hv]);
        }
        auto face = _mesh.add_face(newHfVertices);
        Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
        _mesh.data(face).set_normal_q(n);
        _mesh.property(_cmp, face) = cmpNotSetIdx;
        _mesh.property(_originalHf, face) = hf;
        _mesh.set_normal(face, n.unaryExpr([](mpq_class x) { return x.get_d(); }));
    }

    if (!OpenMesh::IO::write_mesh(_mesh, "debug/_original.obj")) {
        std::cerr << "write error\n";
        exit(1);
    }
}

std::vector<Vector3q> StarDecompositionBoundary::centers() {
    if (!this->_computed) {
        this->run();
    }

    return std::vector<Vector3q>();
}

std::vector<VolumeMesh> StarDecompositionBoundary::components() {
    if (!this->_computed) {
        this->run();
    }

    return std::vector<VolumeMesh>();
}

void StarDecompositionBoundary::run() {
    std::queue<OpenMesh::FaceHandle> queue;
    queue.push(*_mesh.faces_begin());

    _cmpIdx = 0;
    while (!queue.empty()) {
        OpenMesh::FaceHandle h = queue.front();
        queue.pop();
        if (_mesh.property(_cmp, h) != cmpNotSetIdx) {
            continue;
        }

        add_component(h);

        std::deque<OpenMesh::FaceHandle> candidates;
        for (auto he : _mesh.fh_range(h)) {
            auto oFace = _mesh.opposite_face_handle(he);
            if (oFace.is_valid()) {
                candidates.push_back(oFace);
            }
        }
        while (!candidates.empty()) {
            auto nextH = candidates.front();
            candidates.pop_front();
            if (_mesh.property(_cmp, nextH) != cmpNotSetIdx) {
                continue;
            }

            if (!add_hf_to_cmp(_cmpMeshes[_cmpIdx], nextH)) {
                queue.push(nextH);
                continue;
            }

            _mesh.property(_cmp, nextH) = _cmpIdx;
            for (auto he : _mesh.fh_range(nextH)) {
                auto oFace = _mesh.opposite_face_handle(he);
                if (oFace.is_valid() && _mesh.property(_cmp, oFace) == cmpNotSetIdx) {
                    // Check if face can be added between to existing faces
                    bool between = false;
                    for (auto fh : _mesh.fh_range(oFace)) {
                        auto hv = _mesh.from_vertex_handle(fh);
                        auto hvTo = _mesh.to_vertex_handle(fh);
                        if (fh != he && _cmpVertexMap[hv].is_valid() && _cmpVertexMap[hvTo].is_valid()) {
                            auto h = _cmpMeshes[_cmpIdx].find_halfedge(_cmpVertexMap[hv], _cmpVertexMap[hvTo]);
                            if (h.is_valid()) {
                                between = true;
                            }
                            break;
                        }
                    }

                    // Prefer faces that can be added between
                    if (between) {
                        candidates.push_front(oFace);
                    } else {
                        candidates.push_back(oFace);
                    }
                }
            }
        }

        _cmpMeshes[_cmpIdx].garbage_collection();
        {
            char buffer[100];
            std::sprintf(buffer, "debug/cmp_%d.obj", _cmpIdx);
            if (!OpenMesh::IO::write_mesh(_cmpMeshes[_cmpIdx], buffer)) {
                std::cerr << "write error\n";
                exit(1);
            }
        }
        std::cout << "Cmp " << _cmpIdx << " amount faces: " << _cmpMeshes[_cmpIdx].faces().to_vector().size() << std::endl;
    }

    std::cout << "Num components: " << _cmpIdx << std::endl;

    this->_computed = true;
    if (this->_wasVolumeMesh) {
        auto cmp = _originalMesh.property<Cell, int>("cmp", -1);
        for (auto f : _mesh.faces()) {
            auto oF = _originalMesh.opposite_halfface_handle(_mesh.property(_originalHf, f));
            cmp[_originalMesh.incident_cell(oF)] = _mesh.property(_cmp, f);
        }
    }
}

Mesh StarDecompositionBoundary::add_component(OpenMesh::FaceHandle startF) {
    _cmpVertexMap = std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle>();

    Mesh mesh;
    std::vector<OpenMesh::VertexHandle> newHfVertices;
    for (auto hv : _mesh.fv_range(startF)) {
        auto hvQ = _mesh.data(hv).point_q();
        auto newVertex = mesh.add_vertex(Mesh::Point(hvQ[0].get_d(), hvQ[1].get_d(), hvQ[2].get_d()));
        mesh.data(newVertex).set_point_q(hvQ);
        _cmpVertexMap[hv] = newVertex;
        newHfVertices.push_back(newVertex);
    }

    auto face = mesh.add_face(newHfVertices);
    mesh.set_normal(face, _mesh.data(startF).normal_q().unaryExpr([](mpq_class x) { return x.get_d(); }));
    mesh.data(face).set_normal_q(_mesh.data(startF).normal_q());

    // TODO: Check if added faces cut any faces in the original mesh
    auto opposite = get_opposite_face(_mesh, startF);
    if (opposite.first.is_valid()) {
        auto point = opposite.second;
        auto newVertex = mesh.add_vertex(Mesh::Point(point[0].get_d(), point[1].get_d(), point[2].get_d()));
        mesh.data(newVertex).set_point_q(point);
        _cmpFixV = std::make_pair(opposite.first, newVertex);
        for (auto he : mesh.fh_range(face)) {
            auto ohe = mesh.opposite_halfedge_handle(he);
            auto of = mesh.add_face(mesh.from_vertex_handle(ohe), mesh.to_vertex_handle(ohe), newVertex);
            mesh.set_normal(of, mesh.calc_normal(of));
        }
    } else {
        std::cout << "Opposite face not found" << std::endl;
        exit(1);
    }

    _cmpIdx = _cmpMeshes.size();
    _cmpMeshes.push_back(mesh);

    _mesh.property(_cmp, startF) = _cmpIdx;
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
bool StarDecompositionBoundary::add_hf_to_cmp(Mesh& mesh, OpenMesh::FaceHandle& hf) {
    // TODO: Check which edge connects to the component
    // TODO: Use connecting edge to remove an temporary face
    // TODO: Connect open boundaries to the component
    // TODO: Check if opposite faces should be expanded
    OpenMesh::VertexHandle newHfVertex;
    std::vector<OpenMesh::VertexHandle> triangle;
    int deletedHelperFaces = 0;
    OpenMesh::VertexHandle fixedVertex = _cmpFixV.second;
    TxDeleteMesh txMesh(mesh);

    for (auto fh : _mesh.fh_range(hf)) {
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

    for (auto fh : _mesh.fh_range(hf)) {
        auto hv = _mesh.from_vertex_handle(fh);
        if (!_cmpVertexMap[hv].is_valid()) {
            auto hvQ = _mesh.data(hv).point_q();
            auto newVertex = mesh.add_vertex(Mesh::Point(hvQ[0].get_d(), hvQ[1].get_d(), hvQ[2].get_d()));
            mesh.data(newVertex).set_point_q(hvQ);
            _cmpVertexMap[hv] = newVertex;
            newHfVertex = hv;
        }

        triangle.push_back(_cmpVertexMap[hv]);
    }

    bool shouldCheck = false;
    /*if (hf == _cmpFixV.first && newHfVertex.is_valid()) {
        for (auto vf : mesh.vf_range(fixedVertex)) {
            std::vector<OpenMesh::VertexHandle> aTri;
            for (auto v : mesh.fv_range(vf)) {
                if (v == fixedVertex) {
                    aTri.push_back(_cmpVertexMap[newHfVertex]);
                } else {
                    aTri.push_back(v);
                }
            }
            txMesh.delete_face(vf);

            OpenMesh::FaceHandle face = mesh.add_face(aTri);
            mesh.set_normal(face, mesh.calc_normal(face));
            if (mesh.status(fixedVertex).deleted()) {
                break;
            }
        }

        OpenMesh::FaceHandle face = txMesh.add_face(triangle);
        mesh.set_normal(face, _mesh.data(hf).normal_q().unaryExpr([](mpq_class x) { return x.get_d(); }));
    } else */if (deletedHelperFaces != 0 && (deletedHelperFaces != 1 || newHfVertex.is_valid())) {
        OpenMesh::FaceHandle face = txMesh.add_face(triangle);
        mesh.set_normal(face, _mesh.data(hf).normal_q().unaryExpr([](mpq_class x) { return x.get_d(); }));
        for (auto fh : mesh.fh_range(face)) {
            if (mesh.opposite_halfedge_handle(fh).is_boundary()) {
                auto f = txMesh.add_face(mesh.to_vertex_handle(fh), mesh.from_vertex_handle(fh), fixedVertex);
                mesh.set_normal(f, mesh.calc_normal(f));
            }
        }
        shouldCheck = true;
    }

    if (!shouldCheck || !has_valid_center(mesh) || hf == _cmpFixV.first) {
        if (newHfVertex.is_valid()) {
            _cmpVertexMap.erase(newHfVertex);
        }

        txMesh.revert();

        return false;
    }/* else if (hf == _cmpFixV.first) {
        _cmpFixV = std::make_pair(OpenMesh::FaceHandle(), _cmpVertexMap[newHfVertex]);
    }*/


    return true;
}

std::pair<OpenMesh::FaceHandle, Vector3q> StarDecompositionBoundary::get_opposite_face(Mesh& mesh, OpenMesh::FaceHandle& origin) {
    Vector3q vPos = get_face_center(mesh, origin);
    auto n = mesh.data(origin).normal_q();

    double distResult = -1;
    Vector3q pResult;
    OpenMesh::FaceHandle opposite;

    // Find cutting boundary face
    for (auto f : mesh.faces()) {
        if (f == origin) {
            continue;
        }

        auto faceNormal = mesh.data(f).normal_q();
        auto div = faceNormal.transpose() * n;
        if (div[0] == 0) {
            continue;
        }

        std::vector<Vector3q> triangle;
        for (auto hv : mesh.fv_range(f)) {
            triangle.push_back(mesh.data(hv).point_q());
        }

        // Calculate intersection
        auto c = faceNormal.transpose() * mesh.data(mesh.to_vertex_handle(mesh.halfedge_handle(f))).point_q();
        auto nominator = faceNormal.transpose() * vPos - c;
        auto denominator = faceNormal.transpose() * n;
        auto r = (-nominator / denominator)[0];
        if (r > 0) {
            continue;
        }

        auto p = vPos + r * n;
        if (!point_in_triangle(p, triangle[0], triangle[1], triangle[2])) {
            continue;
        }

        double dist = (vPos - p).unaryExpr([](mpq_class x) { return x.get_d(); }).norm();
        if (dist == 0) {
            continue;
        }

        if (distResult == -1 || dist < distResult) {
            distResult = dist;
            pResult = p;
            opposite = f;
        }
    }

    return std::make_pair(opposite, pResult);
}

// https://math.stackexchange.com/questions/51326/determining-if-an-arbitrary-point-lies-inside-a-triangle-defined-by-three-points
// https://gdbooks.gitbooks.io/3dcollisions/content/Chapter4/point_in_triangle.html
bool StarDecompositionBoundary::point_in_triangle(Vector3q p, Vector3q v0, Vector3q v1, Vector3q v2) {
    v0 -= p;
    v1 -= p;
    v2 -= p;

    auto u = v1.cross(v2);
    auto v = v2.cross(v0);
    auto w = v0.cross(v1);

    if (u.dot(v) < 0 || u.dot(w) < 0) {
        return false;
    }

    return true;
}

bool StarDecompositionBoundary::has_valid_center(Mesh& mesh) {
    std::vector<Eigen::Vector3d> positions;
    std::vector<Eigen::Vector3d> normals;
    for (auto face : mesh.faces()) {
        if (face.is_valid() && !face.deleted()) {
            positions.push_back(mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(face))));
            normals.push_back(mesh.normal(face).normalized());
        }
    }

    Vector3q new_center = kernel_chebyshev_center(positions, normals).cast<mpq_class>();
    for (auto face : mesh.faces()) {
        std::vector<Vector3q> triangle;
        for (auto fv : mesh.fv_range(face)) {
            triangle.push_back(mesh.data(fv).point_q());
        }
        Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
        if ((triangle[0] - new_center).dot(n) <= 0) {
            return false;
        }
    }

    return true;
}

Vector3q StarDecompositionBoundary::get_face_center(Mesh& mesh, OpenMesh::FaceHandle& hf) {
    auto v0 = mesh.to_vertex_handle(mesh.halfedge_handle(hf));
    auto v1 = mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(hf)));
    auto v2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.next_halfedge_handle(mesh.halfedge_handle(hf))));

    return (mesh.data(v0).point_q() + mesh.data(v1).point_q() + mesh.data(v2).point_q()) / 3;
}

/*
OpenMesh::FaceHandle StarDecompositionBoundary::check_intersecting(OpenMesh::FaceHandle hf) {
    auto v0 = _mesh.to_vertex_handle(_mesh.halfedge_handle(hf));
    auto v1 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(_mesh.halfedge_handle(hf)));
    auto v2 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(_mesh.next_halfedge_handle(_mesh.halfedge_handle(hf))));

    float p0[3] = { static_cast<float>(_mesh.point(v0)[0]), static_cast<float>(_mesh.point(v0)[1]), static_cast<float>(_mesh.point(v0)[2]) };
    float p1[3] = { static_cast<float>(_mesh.point(v1)[0]), static_cast<float>(_mesh.point(v1)[1]), static_cast<float>(_mesh.point(v1)[2]) };
    float p2[3] = { static_cast<float>(_mesh.point(v2)[0]), static_cast<float>(_mesh.point(v2)[1]), static_cast<float>(_mesh.point(v2)[2]) };
    
    for (auto f : _mesh.faces()) {

    }
    // tri_tri_intersect(_mesh.point(v0).array(), _mesh.point(v1).array(), _mesh.point(v2).array(), _mesh.point(v0).array(), _mesh.point(v1).array(), _mesh.point(v2).array());
}
*/

/*
bool StarDecompositionBoundary::add_halfface(std::set<Halfface>& cmpHf, Halfface& hf) {
    cmpHf.insert(hf);

    if (cmpHf.size() == 1) {
        return true;
    }

    Mesh mesh;
    OpenMesh::FPropHandleT<Vector3q> normal;
    OpenMesh::VPropHandleT<Vector3q> Q;
    mesh.add_property(normal);
    mesh.add_property(Q);
    std::map<Vertex, OpenMesh::SmartVertexHandle> vertexMap;
    for (auto h : cmpHf) {
        std::vector<OpenMesh::SmartVertexHandle> newHfVertices;
        for (auto hv : _mesh.halfface_vertices(h)) {
            if (!vertexMap[hv].is_valid()) {
                auto newVertex = mesh.add_vertex(Mesh::Point(_Q[hv][0].get_d(), _Q[hv][1].get_d(), _Q[hv][2].get_d()));
                vertexMap[hv] = newVertex;
                mesh.property(Q, newVertex) = _Q[hv];
            }
            newHfVertices.push_back(vertexMap[hv]);
        }
        auto face = mesh.add_face(newHfVertices);
        mesh.property(normal, face) = _normal[h];
        mesh.set_normal(face, _normal[h].unaryExpr([](mpq_class x) { return x.get_d(); }));
    }

    // Expand boundary to next opposite halfface
    for (auto heh : mesh.halfedges()) {
        if (!mesh.is_boundary(heh)) {
            continue;
        }

        auto v = mesh.to_vertex_handle(heh);
        auto n = mesh.property(normal, mesh.opposite_face_handle(heh));
        std::vector<std::pair<OpenVolumeMesh::HalfFaceHandle, std::pair<Vector3q, double>>> opposite;
        // Find cutting boundary face
        for (auto hf : _boundary) {
            auto faceNormal = _normal[hf];
            auto div = faceNormal.transpose() * n;
            if (div[0] == 0) {
                continue;
            }

            std::vector<Eigen::Vector3d> triangle;
            for (auto hv : _mesh.halfface_vertices(hf)) {
                triangle.push_back(_Q[hv].unaryExpr([](mpq_class x) { return x.get_d(); }));
            }
            auto vPos = mesh.property(Q, v);

            // Calculate intersection
            auto c = faceNormal.transpose() * _Q[*_mesh.halfface_vertices(hf).first];
            auto nominator = faceNormal.transpose() * vPos - c;
            auto denominator = faceNormal.transpose() * n;
            double r = (-nominator / denominator)[0].get_d();
            if (r < 0) {
                continue;
            }

            auto p = vPos + r * n;
            Eigen::Vector3d pD = p.unaryExpr([](mpq_class x) { return x.get_d(); });
            // Calculate baricentric coordinates
            auto A = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]).norm() / 2;
            auto Aa = (triangle[0] - pD).cross(triangle[1] - pD).norm() / 2;
            auto Ab = (triangle[1] - pD).cross(triangle[2] - pD).norm() / 2;
            auto Ac = (triangle[0] - pD).cross(triangle[2] - pD).norm() / 2;

            if (Aa / A + Ab / A + Ac / A != 1) {
                continue;
            }

            double dist = (vPos - p).unaryExpr([](mpq_class x) { return x.get_d(); }).norm();
            if (dist == 0) {
                continue;
            }

            auto result = std::make_pair(hf, std::make_pair(p, dist));
            if (opposite.size() == 0) {
                opposite.push_back(result);
                continue;
            }

            if (dist < opposite[0].second.second) {
                opposite[0] = result;
            }
        }

        if (opposite.size() > 0) {
            // std::cout << opposite[0].second.second << std::endl;
            auto newVPos = opposite[0].second.first;
            auto newV = mesh.add_vertex(Mesh::Point(newVPos[0].get_d(), newVPos[1].get_d(), newVPos[2].get_d()));
            mesh.property(Q, newV) = newVPos;
            auto face = mesh.add_face(mesh.from_vertex_handle(heh), v, newV);
            mesh.property(normal, face) = mesh.calc_normal(face).cast<mpq_class>();
            mesh.set_normal(face, mesh.calc_normal(face));
        }
    }

    {
        char buffer[100];
        std::sprintf(buffer, "debug/cmp_%d_before_close.obj", _cmpIdx);
        if (!OpenMesh::IO::write_mesh(mesh, buffer)) {
            std::cerr << "write error\n";
            exit(1);
        }
    }

    // Close boundary
    bool gotToEnd = false;
    int addedFaces = 0;
    while (!gotToEnd) {
        gotToEnd = true;
        for (auto startHeh : mesh.halfedges()) {
            if (!mesh.is_boundary(startHeh)) {
                continue;
            }

            auto heh = mesh.next_halfedge_handle(startHeh);
            if (heh == startHeh) {
                std::cout << "Error: startHeh == heh" << std::endl;
                continue;
            }

            if (!mesh.is_boundary(heh)) {
                std::cout << "Error: heh non boundary" << std::endl;
                continue;
            }

            auto newEdge = mesh.find_halfedge(mesh.to_vertex_handle(heh), mesh.from_vertex_handle(startHeh));
            if (newEdge.is_valid() && !newEdge.is_boundary()) {
                newEdge = mesh.find_halfedge(mesh.from_vertex_handle(startHeh), mesh.to_vertex_handle(heh));
                if (newEdge.is_valid() && !newEdge.is_boundary()) {
                    continue;
                }
            }

            auto ndA = mesh.property(normal, mesh.opposite_face_handle(startHeh)).unaryExpr([](mpq_class x) { return x.get_d(); });
            auto ndB = mesh.property(normal, mesh.opposite_face_handle(heh)).unaryExpr([](mpq_class x) { return x.get_d(); });

            std::vector<OpenMesh::VertexHandle> triangle;
            triangle.push_back(mesh.from_vertex_handle(startHeh));
            triangle.push_back(mesh.to_vertex_handle(startHeh));
            triangle.push_back(mesh.to_vertex_handle(heh));

            // TODO: Do not add face if it cuts through face in _mesh

            Vector3q n = (mesh.property(Q, triangle[1]) - mesh.property(Q, triangle[0])).cross(mesh.property(Q, triangle[2]) - mesh.property(Q, triangle[0]));
            Eigen::Vector3d nd = n.unaryExpr([](mpq_class x) { return x.get_d(); });
            if (nd.norm() == 0) {
                continue;
            }

            double s1 = nd.normalized().dot(mesh.normal(mesh.opposite_face_handle(startHeh)).normalized());
            double s2 = nd.normalized().dot(mesh.normal(mesh.opposite_face_handle(heh)).normalized());
            if (s1 > 0 && s2 > 0) {
                continue;
            }
            
            auto face = mesh.add_face(triangle);
            if (!face.is_valid()) {
                throw std::runtime_error("Error: could not add face to component mesh");
                continue;
            }

            mesh.property(normal, face) = n;
            mesh.set_normal(face, mesh.calc_normal(face));
            gotToEnd = false;
            addedFaces++;

            {
                char buffer[100];
                std::sprintf(buffer, "debug/cmp_%d_close_%d.obj", _cmpIdx, addedFaces);
                if (!OpenMesh::IO::write_mesh(mesh, buffer)) {
                    std::cerr << "write error\n";
                    exit(1);
                }
            }

            break;
        }
    }

    bool new_center_valid = has_valid_center(mesh);
    if (!new_center_valid) {
        char buffer[100];
        std::sprintf(buffer, "debug/cmp_%d_invalid.obj", _cmpIdx);
        if (!OpenMesh::IO::write_mesh(mesh, buffer)) {
            std::cerr << "write error\n";
            exit(1);
        }

        cmpHf.erase(hf);
    } else {
        char buffer[100];
        std::sprintf(buffer, "debug/cmp_%d_valid.obj", _cmpIdx);
        if (!OpenMesh::IO::write_mesh(mesh, buffer)) {
            std::cerr << "write error\n";
            exit(1);
        }
    }

    return new_center_valid;
}
*/
