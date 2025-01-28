#include "sd_boundary.h"
#include "lp.h"
#include "tritri.h"
#include "assertion.h"
#include <Eigen/src/Core/Matrix.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <deque>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <queue>

// #define FACE_SELECTION_ORDERED 1

const int cmpNotSetIdx = -1;

StarDecompositionBoundary::StarDecompositionBoundary(Mesh& m) : _computed(false), _wasVolumeMesh(false), _mesh(m) {
    _mesh.add_property(_cmp);
    _mesh.add_property(_originalHf);
    for (auto v : _mesh.vertices()) {
        _mesh.data(v).set_point_q(_mesh.point(v).cast<mpq_class>());
    }

    for (auto f : _mesh.faces()) {
        _mesh.update_normal_q(f);
        _mesh.property(_cmp, f) = cmpNotSetIdx;
    }

    if (!OpenMesh::IO::write_mesh(_mesh, "debug/original.obj")) {
        std::cerr << "write error\n";
        exit(1);
    }
}

StarDecompositionBoundary::StarDecompositionBoundary(VolumeMesh& mesh) : _computed(false), _wasVolumeMesh(true), _originalMesh(mesh) {
    _mesh.add_property(_cmp);
    _mesh.add_property(_originalHf);

    auto Q = mesh.property<Vertex, Vector3q>("Q");
    std::map<Vertex, OpenMesh::VertexHandle> vMap;
    for (auto v : mesh.vertices()) {
        vMap[v] = _mesh.add_vertex_q(Q[v]);
    }

    for (auto hf : mesh.boundary_halffaces()) {
        std::vector<OpenMesh::VertexHandle> newHfVertices;
        for (auto hv : mesh.halfface_vertices(hf)) {
            newHfVertices.push_back(vMap[hv]);
        }
        auto face = _mesh.add_face(newHfVertices);
        _mesh.update_normal_q(face);
        _mesh.property(_cmp, face) = cmpNotSetIdx;
        _mesh.property(_originalHf, face) = hf;
    }

    if (!OpenMesh::IO::write_mesh(_mesh, "debug/original.obj")) {
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

struct QueueCandidate
{
    double priority{};
    OpenMesh::FaceHandle data{};

    friend bool operator>(QueueCandidate const& lhs, QueueCandidate const& rhs)
    {
        return lhs.priority > rhs.priority;
    }
};

void StarDecompositionBoundary::run() {
#ifndef DEBUG
    srand(0);
#endif

    _cmpIdx = 0;
    int startOffset = 0;
    int resets = 0;
    int bestSinceReset = 0;
    while (_mesh.faces_empty() == false) {
#ifdef FACE_SELECTION_ORDERED
        auto nextFace = _mesh.faces_begin() + startOffset;
        if (nextFace == _mesh.faces_end()) {
            std::cout << "Decomposition failed" << std::endl;
            break;
        }
        OpenMesh::FaceHandle h = *nextFace;
        if (_mesh.property(_cmp, h) != cmpNotSetIdx) {
            continue;
        }
#else
        // Add random face
        OpenMesh::FaceHandle h = *(_mesh.faces_begin() + (rand() % _mesh.n_faces()));
#endif

        add_component(h);
        if (has_valid_center(_mesh).first) {
            for (auto f : _mesh.faces()) {
                _mesh.property(_cmp, f) = cmpNotSetIdx;
                _cmpMeshes[_cmpIdx] = _mesh;
            }
            break;
        }

        std::priority_queue<QueueCandidate, std::vector<QueueCandidate>, std::greater<QueueCandidate>> candidates;
        std::set<OpenMesh::FaceHandle> visited;
        for (auto he : _mesh.fh_range(h)) {
            auto oFace = _mesh.opposite_face_handle(he);
            if (oFace.is_valid()) {
                candidates.push({ 0, oFace});
            }
        }

        int allowMove = 0;
        Vector3q beforeMove = _currentCmp->data(_cmpFixV.second).point_q();
        while (!candidates.empty()) {
            OpenMesh::FaceHandle nextH;
            nextH = candidates.top().data;
            candidates.pop();
            bool alreadyChecked = _mesh.property(_cmp, nextH) != cmpNotSetIdx || visited.find(nextH) != visited.end();
            if (alreadyChecked || !add_face_to_cmp(*_currentCmp, nextH)) {
                visited.insert(nextH);
                if (candidates.empty() && allowMove > 0) {
                    beforeMove = _currentCmp->data(_cmpFixV.second).point_q();
                    if (!move_fix_vertex(*_currentCmp)) {
                        break;
                    }

                    for (auto f : _currentCmp->vf_range(_cmpFixV.second)) {
                        auto he = _currentCmp->halfedge_handle(f);
                        for (auto fhe : _currentCmp->fh_range(f)) {
                            auto fromV = _currentCmp->from_vertex_handle(fhe);
                            auto toV = _currentCmp->to_vertex_handle(fhe);
                            if (fromV == _cmpFixV.second || toV == _cmpFixV.second) {
                                continue;
                            }

                            auto mH = _mesh.find_halfedge(_meshVertexMap[fromV], _meshVertexMap[toV]);
                            auto dist = _mesh.face_center(_mesh.face_handle(mH)) - _cmpCenter;
                            double distAbs = std::abs(dist[0].get_d()) + std::abs(dist[1].get_d()) + std::abs(dist[2].get_d());
                            candidates.push({ distAbs, _mesh.face_handle(mH) });
                        }
                    }
                    visited.clear();
                    allowMove = 0;
                } else if (candidates.empty() && allowMove == 0 && beforeMove != _currentCmp->data(_cmpFixV.second).point_q()) {
                    _currentCmp->move(_cmpFixV.second, beforeMove);
                }
                continue;
            }

            visited.clear();
            allowMove++;
            _mesh.property(_cmp, nextH) = _cmpIdx;
            for (auto he : _mesh.fh_range(nextH)) {
                auto oFace = _mesh.opposite_face_handle(he);
                if (oFace.is_valid() && _mesh.property(_cmp, oFace) == cmpNotSetIdx) {
                    auto dist = _mesh.face_center(oFace) - _cmpCenter;
                    double distAbs = std::abs(dist[0].get_d()) + std::abs(dist[1].get_d()) + std::abs(dist[2].get_d());
                    candidates.push({ distAbs, oFace });
                }
            }
        }

        int numCmpFacesTotal = _currentCmp->faces().to_vector().size();
        int numCmpFaces = 0;
        for (auto f : _mesh.faces()) {
            if (_mesh.property(_cmp, f) == _cmpIdx) {
                numCmpFaces++;
            }
        }
        int amount = (-numCmpFaces + (numCmpFacesTotal - numCmpFaces));
        if ((resets > 5 && amount <= 0.3 * bestSinceReset)) {
            std::cout << "Cmp " << _cmpIdx << " amount faces: " << _currentCmp->faces().to_vector().size() << " full mesh size change: " << (-numCmpFaces + (numCmpFacesTotal - numCmpFaces));
            std::cout << " resets: " << resets << " full mesh size: " << _mesh.n_faces() << std::endl;
            apply_current_component();
            {
                char buffer[100];
                std::sprintf(buffer, "debug/original_%d.obj", _cmpIdx);
                if (!OpenMesh::IO::write_mesh(_mesh, buffer)) {
                    std::cerr << "write error\n";
                    exit(1);
                }
            }
            startOffset = 0;
            resets = 0;
            bestSinceReset = 0;
        } else {
            for (auto f : _mesh.faces()) {
                if (_mesh.property(_cmp, f) == _cmpIdx) {
                    _mesh.property(_cmp, f) = cmpNotSetIdx;
                }
            }
        
            _cmpMeshes.pop_back();
            startOffset++;
            resets++;

            if (amount < bestSinceReset) {
                bestSinceReset = (-numCmpFaces + (numCmpFacesTotal - numCmpFaces));
            }
        }

        _currentCmp->garbage_collection();
        {
            char buffer[100];
            std::sprintf(buffer, "debug/cmp_%d.obj", _cmpIdx);
            if (!OpenMesh::IO::write_mesh(*_currentCmp, buffer)) {
                std::cerr << "write error\n";
                exit(1);
            }
        }
    }

    std::cout << "Num components: " << _cmpIdx + 1 << std::endl;
    {
        char buffer[100];
        std::sprintf(buffer, "debug/cmp_%d.obj", _cmpIdx);
        if (!OpenMesh::IO::write_mesh(*_currentCmp, buffer)) {
            std::cerr << "write error\n";
            exit(1);
        }
    }

    this->_computed = true;
    if (this->_wasVolumeMesh) {
        auto cmp = _originalMesh.property<Cell, int>("cmp", -1);
        for (auto f : _mesh.faces()) {
            auto oF = _originalMesh.opposite_halfface_handle(_mesh.property(_originalHf, f));
            cmp[_originalMesh.incident_cell(oF)] = _mesh.property(_cmp, f);
        }
    }
}

void StarDecompositionBoundary::apply_current_component() {
    Mesh& cmpMesh = *_currentCmp;
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
            std::cout << "Improved fix vertex pos" << std::endl;
            break;
        }
    }

    // for (auto f : cmpMesh.vf_range(_cmpFixV.second).to_vector()) {
        // std::cout << "Fix vertex face: " << f.idx() << std::endl;
        // cmpMesh.split(f);
    // }

    double maxArea = -1;
    auto fixVertex = _mesh.add_vertex_q(cmpMesh.data(_cmpFixV.second).point_q());
    for (auto f : _mesh.faces()) {
        if (_mesh.property(_cmp, f) == _cmpIdx) {
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
            _mesh.property(_cmp, f) = cmpNotSetIdx;
        }
    }

    std::function<void(OpenMesh::VertexHandle)> recShrinkFaces;
    recShrinkFaces = [this, maxArea, &recShrinkFaces](OpenMesh::VertexHandle v) {
        for (auto f : _mesh.vf_range(v)) {
            if (_mesh.calc_face_area(f) > maxArea) {
                std::cout << "Shrink face" << std::endl;
                auto center = _mesh.face_center(f);
                auto v = _mesh.split(f, center.unaryExpr([](mpq_class x) { return x.get_d(); }));
                _mesh.data(v).set_point_q(center);
                for (auto vf : _mesh.vf_range(v)) {
                    _mesh.update_normal_q(vf);
                }
                recShrinkFaces(v);
            }
        }
    };
    recShrinkFaces(fixVertex);

    _mesh.delete_isolated_vertices();
    _mesh.garbage_collection();
}

bool StarDecompositionBoundary::move_fix_vertex(Mesh& mesh) {
    auto opposite = get_opposite_face(mesh, _cmpFixV.second);
    if (!opposite.is_valid()) {
        std::cout << "Opposite face not found" << std::endl;
        return false;
    }

    Vector3q vPos = mesh.data(_cmpFixV.second).point_q();
    auto n = -mesh.calc_normal(_cmpFixV.second).cast<mpq_class>();
    auto faceNormal = _mesh.data(opposite).normal_q();

    // Calculate intersection
    auto c = faceNormal.transpose() * _mesh.data(_mesh.to_vertex_handle(_mesh.halfedge_handle(opposite))).point_q();
    auto nominator = faceNormal.transpose() * vPos - c;
    auto denominator = faceNormal.transpose() * n;
    auto r = (-nominator / denominator)[0];
    r /= 2;
    Vector3q p = vPos + r * n;

    p = p.unaryExpr([](mpq_class x) { return x.get_d(); }).cast<mpq_class>();
    for (int i = 0; i < 3; i++) {
        if (move_vertex_to(mesh, _cmpFixV.second, p)) {
            _cmpFixV.first = opposite;
            return true;
        } else {
            r /= 2;
            p = vPos + r * n;
            p = p.unaryExpr([](mpq_class x) { return x.get_d(); }).cast<mpq_class>();
        }
    }

    return false;
}

bool StarDecompositionBoundary::move_vertex_to(Mesh& mesh, OpenMesh::VertexHandle& moveV, const Vector3q& p) {
    for (auto vf : mesh.vf_range(moveV)) {
        std::vector<OpenMesh::VertexHandle> vertices;
        for (auto v : mesh.fv_range(vf)) {
            if (moveV != v) {
                vertices.push_back(v);
            }
        }
        std::vector<Vector3q> t = { mesh.data(vertices[0]).point_q(), mesh.data(vertices[1]).point_q(), p };
        if (triangles_intersect(t, { _meshVertexMap[vertices[0]], _meshVertexMap[vertices[1]] })) {
            return false;
        }
    }

    Vector3q prevPos = mesh.data(moveV).point_q();
    mesh.move(_cmpFixV.second, p);
    if (!has_valid_center(mesh).first) {
        mesh.move(_cmpFixV.second, prevPos);
        return false;
    }

    return true;
}

Mesh StarDecompositionBoundary::add_component(OpenMesh::FaceHandle startF) {
    _cmpVertexMap = std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle>();
    _meshVertexMap = std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle>();

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

    _cmpIdx = _cmpMeshes.size();
    _cmpMeshes.push_back(mesh);
    _currentCmp = &_cmpMeshes[_cmpIdx];

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
bool StarDecompositionBoundary::add_face_to_cmp(Mesh& mesh, OpenMesh::FaceHandle& newFace) {
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
    bool illegalTriangle = false;
    if (deletedHelperFaces != 0 && (deletedHelperFaces != 1 || newFaceVertex.is_valid())) {
        OpenMesh::FaceHandle face = txMesh.add_face(triangle);
        mesh.set_normal_q(face, _mesh.data(newFace).normal_q());
        for (auto fh : mesh.fh_range(face)) {
            if (mesh.opposite_halfedge_handle(fh).is_boundary()) {
                auto v0 = mesh.to_vertex_handle(fh);
                auto v1 = mesh.from_vertex_handle(fh);
                auto v2 = fixedVertex;
                std::vector<Vector3q> t = { mesh.data(v0).point_q(), mesh.data(v1).point_q(), mesh.data(v2).point_q() };
                if (triangles_intersect(t, { _meshVertexMap[v0], _meshVertexMap[v1] })) {
                    illegalTriangle = true;
                    break;
                }

                auto f = txMesh.add_face(v0, v1, v2);
                mesh.update_normal_q(f);
            }
        }
        shouldCheck = true;
    }

    std::pair<bool, Vector3q> cmpCenter;
    if (shouldCheck && !illegalTriangle) {
        cmpCenter = has_valid_center(mesh);
    }

    if (illegalTriangle || !shouldCheck || !cmpCenter.first) {
        if (newFaceVertex.is_valid()) {
            _meshVertexMap.erase(_cmpVertexMap[newFaceVertex]);
            _cmpVertexMap.erase(newFaceVertex);
        }

        txMesh.revert();

        return false;
    }


    _cmpCenter = cmpCenter.second;
    return true;
}

OpenMesh::FaceHandle StarDecompositionBoundary::get_opposite_face(Mesh& mesh, OpenMesh::FaceHandle& origin) {
    Vector3q vPos = mesh.face_center(origin);
    Vector3q n = mesh.data(origin).normal_q();

    return get_opposite_face(vPos, n);
}

OpenMesh::FaceHandle StarDecompositionBoundary::get_opposite_face(Mesh& mesh, OpenMesh::VertexHandle& origin) {
    Vector3q vPos = mesh.data(origin).point_q();
    auto n = -mesh.calc_normal(origin);

    return get_opposite_face(vPos, n.cast<mpq_class>());
}

OpenMesh::FaceHandle StarDecompositionBoundary::get_opposite_face(Vector3q vPos, Vector3q n) {
    double distResult = -1;
    Vector3q pResult;
    OpenMesh::FaceHandle opposite;

    // Find cutting boundary face
    for (auto f : _mesh.faces()) {
        auto faceNormal = _mesh.data(f).normal_q();
        auto div = faceNormal.transpose() * n;
        if (div[0] == 0) {
            continue;
        }

        // Calculate intersection
        auto c = faceNormal.transpose() * _mesh.data(_mesh.to_vertex_handle(_mesh.halfedge_handle(f))).point_q();
        auto nominator = faceNormal.transpose() * vPos - c;
        auto denominator = faceNormal.transpose() * n;
        auto r = (-nominator / denominator)[0];
        if (r > 0) {
            continue;
        }

        auto p = vPos + r * n;
        if (!_mesh.point_on_face(f, p)) {
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

    return opposite;
}

std::pair<OpenMesh::FaceHandle, Vector3q> StarDecompositionBoundary::get_fix_vertex_pos(Mesh& mesh, OpenMesh::FaceHandle& hf) {
    auto opposite = get_opposite_face(mesh, hf);
    if (!opposite.is_valid()) {
        return std::make_pair(opposite, Vector3q::Zero());
    }

    Vector3q vPos = mesh.face_center(hf);
    auto n = mesh.data(hf).normal_q();
    auto faceNormal = mesh.data(opposite).normal_q();

    // Calculate intersection
    auto c = faceNormal.transpose() * mesh.data(mesh.to_vertex_handle(mesh.halfedge_handle(opposite))).point_q();
    auto nominator = faceNormal.transpose() * vPos - c;
    auto denominator = faceNormal.transpose() * n;
    auto r = (-nominator / denominator)[0];
    r = r / 2;
    Vector3q p = vPos + r * n;

    bool intersects;
    do {
        intersects = false;
        for (auto fh : mesh.fh_range(hf)) {
            auto v0 = mesh.to_vertex_handle(fh);
            auto v1 = mesh.from_vertex_handle(fh);
            std::vector<Vector3q> t = { mesh.data(v0).point_q(), mesh.data(v1).point_q(), p };
            if (triangles_intersect(t, { v0, v1 })) {
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

bool StarDecompositionBoundary::triangles_intersect(std::vector<Vector3q> t, std::vector<OpenMesh::VertexHandle> borderVertices) {
    Vector3q n = (t[1] - t[0]).cross(t[2] - t[0]);
    for (auto face : _mesh.faces()) {
        auto v0 = _mesh.to_vertex_handle(_mesh.halfedge_handle(face));
        auto v1 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(_mesh.halfedge_handle(face)));
        auto v2 = _mesh.to_vertex_handle(_mesh.next_halfedge_handle(_mesh.next_halfedge_handle(_mesh.halfedge_handle(face))));

        // Move connecting vertices slightly to avoid intersection
        bool shareV0 = false, shareV1 = false, shareV2 = false;
        for (auto v : borderVertices) {
            if (v == v0) {
                shareV0 = true;
            } else if (v == v1) {
                shareV1 = true;
            } else if (v == v2) {
                shareV2 = true;
            }
        }
        short sharedVertices = shareV0 + shareV1 + shareV2;
        ASSERT(sharedVertices < 3);

        Vector3q fNormal = _mesh.data(face).normal_q();
        Vector3q v0q = _mesh.data(v0).point_q();
        Vector3q v1q = _mesh.data(v1).point_q();
        Vector3q v2q = _mesh.data(v2).point_q();
        if (shareV0 && shareV1) {
            v0q = v0q + (v2q - v0q) * 1e-6;
            v1q = v1q + (v2q - v1q) * 1e-6;
        } else if (shareV0 && shareV2) {
            v0q = v0q + (v1q - v0q) * 1e-6;
            v2q = v2q + (v1q - v2q) * 1e-6;
        } else if (shareV1 && shareV2) {
            v1q = v1q + (v0q - v1q) * 1e-6;
            v2q = v2q + (v0q - v2q) * 1e-6;
        } else if (shareV0) {
            v0q = v0q + (v2q - v0q) * 1e-6 + (v1q - v0q) * 1e-6;
        } else if (shareV1) {
            v1q = v1q + (v2q - v1q) * 1e-6 + (v0q - v1q) * 1e-6;
        } else if (shareV2) {
            v2q = v2q + (v1q - v2q) * 1e-6 + (v0q - v2q) * 1e-6;
        }

        if (tri_tri_intersect(t[0], t[1], t[2], n, v0q, v1q, v2q, fNormal)) {
            return true;
        }
    }

    return false;
}

std::pair<bool, Vector3q> StarDecompositionBoundary::has_valid_center(Mesh& mesh) {
    std::vector<Eigen::Vector3d> positions;
    std::vector<Eigen::Vector3d> normals;
    for (auto face : mesh.faces()) {
        if (face.is_valid() && !face.deleted()) {
            positions.push_back(mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(face))));
            normals.push_back(mesh.normal(face).normalized());
        }
    }

    Vector3q newCenter = kernel_chebyshev_center(positions, normals).cast<mpq_class>();
    for (auto face : mesh.faces()) {
        std::vector<Vector3q> triangle;
        for (auto fv : mesh.fv_range(face)) {
            triangle.push_back(mesh.data(fv).point_q());
        }
        Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
        if ((triangle[0] - newCenter).dot(n) <= 0) {
            return std::make_pair(false, Vector3q::Zero());
        }
    }

    return std::make_pair(true, newCenter);
}
