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
#include <queue>

#define FACE_SELECTION_ORDERED 1

const int cmpNotSetIdx = -1;

StarDecompositionBoundaryLp::StarDecompositionBoundaryLp(Mesh& m) : _computed(false), _mesh(m) {
    _mesh.add_property(_cmp);
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

StarDecompositionBoundaryLp::StarDecompositionBoundaryLp(VolumeMesh& mesh) : _computed(false), _originalMesh(mesh) {
    _mesh.add_property(_cmp);

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
    }

    if (!OpenMesh::IO::write_mesh(_mesh, "debug/original.obj")) {
        std::cerr << "write error\n";
        exit(1);
    }
}

std::vector<Vector3q> StarDecompositionBoundaryLp::centers() {
    if (!this->_computed) {
        this->run();
    }

    return std::vector<Vector3q>();
}

std::vector<VolumeMesh> StarDecompositionBoundaryLp::components() {
    if (!this->_computed) {
        this->run();
    }

    std::vector<VolumeMesh> components;
    for (auto mesh : _cmpMeshes) {
        VolumeMesh m;
        if (!retetrahedrize(mesh, m)) {
            std::cerr << "Could not tetrahederize component" << std::endl;
            exit(1);
        }

        components.push_back(m);
    }

    return components;
}

void StarDecompositionBoundaryLp::run() {
#ifndef DEBUG
    srand(0);
#endif

    _cmpIdx = 0;
    int startOffset = 0;
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
        if (_mesh.star_center().first) {
            for (auto f : _mesh.faces()) {
                _mesh.property(_cmp, f) = cmpNotSetIdx;
            }
            _cmpMeshes[_cmpIdx] = _mesh;
            break;
        }

        std::queue<OpenMesh::FaceHandle> candidates;
        std::queue<OpenMesh::FaceHandle> candidates2;
        std::set<OpenMesh::FaceHandle> visited;
        for (auto he : _mesh.fh_range(h)) {
            auto oFace = _mesh.opposite_face_handle(he);
            if (oFace.is_valid()) {
                candidates.push(oFace);
            }
        }

        int allowMove = 0;
        while (!candidates.empty() || !candidates2.empty()) {
            OpenMesh::FaceHandle nextH;
            if (!candidates2.empty()) {
                nextH = candidates2.front();
                candidates2.pop();
            } else {
                nextH = candidates.front();
                candidates.pop();
            }
            bool alreadyChecked = _mesh.property(_cmp, nextH) != cmpNotSetIdx || visited.find(nextH) != visited.end();
            if (alreadyChecked || !add_face_to_cmp(*_currentCmp, nextH)) {
                visited.insert(nextH);
                continue;
            }

            visited.clear();
            allowMove++;
            _mesh.property(_cmp, nextH) = _cmpIdx;
            for (auto he : _mesh.fh_range(nextH)) {
                auto oFace = _mesh.opposite_face_handle(he);
                if (oFace.is_valid() && _mesh.property(_cmp, oFace) == cmpNotSetIdx) {
                    // Check if face can be added between two existing faces
                    bool between = false;
                    for (auto fh : _mesh.fh_range(oFace)) {
                        auto hv = _mesh.from_vertex_handle(fh);
                        auto hvTo = _mesh.to_vertex_handle(fh);
                        if (fh != he && _cmpVertexMap[hv].is_valid() && _cmpVertexMap[hvTo].is_valid()) {
                            auto h = _currentCmp->find_halfedge(_cmpVertexMap[hv], _cmpVertexMap[hvTo]);
                            if (h.is_valid()) {
                                between = true;
                            }
                            break;
                        }
                    }

                    if (between) {
                        candidates2.push(oFace);
                    } else {
                        candidates.push(oFace);
                    }
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

        std::cout << "Cmp " << _cmpIdx << " amount faces: " << _currentCmp->faces().to_vector().size() << " full mesh size: " << (-numCmpFaces + (numCmpFacesTotal - numCmpFaces)) << std::endl;
        apply_current_component();
        {
            char buffer[100];
            std::sprintf(buffer, "debug/original_%d.obj", _cmpIdx);
            if (!OpenMesh::IO::write_mesh(_mesh, buffer)) {
                std::cerr << "write error\n";
                exit(1);
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
        exit(1);
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
}

void StarDecompositionBoundaryLp::apply_current_component() {
    for (auto f : _mesh.faces()) {
        if (_mesh.property(_cmp, f) == _cmpIdx) {
            _mesh.delete_face(f, true);
        }
    }

    Mesh& cmpMesh = *_currentCmp;
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
            _mesh.property(_cmp, f) = cmpNotSetIdx;
        }
    }

    _mesh.delete_isolated_vertices();
    _mesh.garbage_collection();
}

Mesh StarDecompositionBoundaryLp::add_component(const OpenMesh::FaceHandle& startF) {
    _cmpVertexMap = std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle>();
    _meshVertexMap = std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle>();

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

    _cmpIdx = _cmpMeshes.size();
    _cmpMeshes.push_back(mesh);
    _currentCmp = &_cmpMeshes[_cmpIdx];

    _cmpNormal = _mesh.data(startF).normal_q();
    _cmpCenter = _mesh.face_center(startF);
    _cmpFixV = _cmpCenter + _mesh.data(startF).normal_q() * 1e-6;

    _mesh.property(_cmp, startF) = _cmpIdx;
    return mesh;
}

/**
 * Adds a face to an existing triangle mesh
 */
bool StarDecompositionBoundaryLp::add_face_to_cmp(Mesh& mesh, OpenMesh::FaceHandle& newFace) {
    OpenMesh::VertexHandle newFaceVertex;
    std::vector<OpenMesh::VertexHandle> triangle;
    TxDeleteMesh txMesh(mesh);

    int existingEdges = 0;
    int borderVertex = 0;
    for (auto fh : _mesh.fh_range(newFace)) {
        auto hv = _mesh.from_vertex_handle(fh);
        auto hvTo = _mesh.to_vertex_handle(fh);
        if (!_cmpVertexMap[hv].is_valid()) {
            auto newVertex = mesh.add_vertex_q(_mesh.data(hv).point_q());
            _cmpVertexMap[hv] = newVertex;
            _meshVertexMap[newVertex] = hv;
            newFaceVertex = hv;
            borderVertex = triangle.size();
        } else if (_currentCmp->find_halfedge(_cmpVertexMap[hvTo], _cmpVertexMap[hv]).is_valid()) {
            existingEdges++;
        } else {
            borderVertex = triangle.size();
        }

        triangle.push_back(_cmpVertexMap[hv]);
    }

    if (existingEdges == 1 && !newFaceVertex.is_valid()) {
        return false;
    }

    bool valid = true;
    OpenMesh::FaceHandle face = txMesh.add_face(triangle);
    mesh.set_normal_q(face, _mesh.data(newFace).normal_q());
    Vector3q meshNormal = (_cmpNormal * (mesh.n_faces() - 1) / mesh.n_faces()) + (_mesh.data(newFace).normal_q() * (1 / mesh.n_faces()));
    if ((_mesh.data(*_mesh.fv_begin(newFace)).point_q() - _cmpCenter).dot(_mesh.data(newFace).normal_q()) <= 0) {
        auto cmpCenter = has_valid_center(mesh, meshNormal);
        if (cmpCenter.first == INVALID) {
            valid = false;
        } else {
            // TODO: Center needs to be moved
            for (auto h : mesh.halfedges()) {
                if (h.is_boundary()) {
                    auto v0 = mesh.to_vertex_handle(h);
                    auto v1 = mesh.from_vertex_handle(h);
                    std::vector<Vector3q> t = { mesh.data(v0).point_q(), mesh.data(v1).point_q(), cmpCenter.second };
                    if (_mesh.triangle_intersects(t, { _meshVertexMap[v0], _meshVertexMap[v1] }).is_valid()) {
                        valid = false;
                        break;
                    }
                }
            }

            if (valid) {
                _cmpFixV = cmpCenter.second;
            }
        }
    } else {
        for (auto fh : mesh.fh_range(face)) {
            if (mesh.opposite_halfedge_handle(fh).is_boundary()) {
                auto v0 = mesh.to_vertex_handle(fh);
                auto v1 = mesh.from_vertex_handle(fh);
                std::vector<Vector3q> t = { mesh.data(v0).point_q(), mesh.data(v1).point_q(), _cmpFixV };
                if (_mesh.triangle_intersects(t, { _meshVertexMap[v0], _meshVertexMap[v1] }).is_valid()) {
                    valid = false;
                    break;
                }
            }
        }
    }

    if (valid == false) {
        if (newFaceVertex.is_valid()) {
            _meshVertexMap.erase(_cmpVertexMap[newFaceVertex]);
            _cmpVertexMap.erase(newFaceVertex);
        }

        txMesh.revert();
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

    auto result = star_center_close_to(normal.unaryExpr([](mpq_class x) { return x.get_d(); }), positions, normals);
    if (result.first == INVALID) {
        return std::make_pair(result.first, Vector3q::Zero());
    }

    Vector3q newCenter = result.second.cast<mpq_class>();
    return std::make_pair(result.first, newCenter);
}
