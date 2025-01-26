#include "sd_boundary_lp.h"
#include "OpenMesh/Core/Mesh/Handles.hh"
#include "lp.h"
#include "tritri.h"
#include "assertion.h"
#include <Eigen/src/Core/Matrix.h>
#include <cmath>
#include <cstdio>
#include <map>
#include <deque>
#include <OpenMesh/Core/IO/MeshIO.hh>

#define FACE_SELECTION_ORDERED 1

const int cmpNotSetIdx = -1;

StarDecompositionBoundaryLp::StarDecompositionBoundaryLp(Mesh& m) : _computed(false), _wasVolumeMesh(false), _mesh(m) {
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

StarDecompositionBoundaryLp::StarDecompositionBoundaryLp(VolumeMesh& mesh) : _computed(false), _wasVolumeMesh(true), _originalMesh(mesh) {
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

    return std::vector<VolumeMesh>();
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
        if (has_valid_center(_mesh).first) {
            for (auto f : _mesh.faces()) {
                _mesh.property(_cmp, f) = cmpNotSetIdx;
                _cmpMeshes[_cmpIdx] = _mesh;
            }
            break;
        }

        std::deque<OpenMesh::FaceHandle> candidates;
        std::set<OpenMesh::FaceHandle> visited;
        for (auto he : _mesh.fh_range(h)) {
            auto oFace = _mesh.opposite_face_handle(he);
            if (oFace.is_valid()) {
                candidates.push_back(oFace);
            }
        }

        int allowMove = 0;
        while (!candidates.empty()) {
            auto nextH = candidates.front();
            candidates.pop_front();
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
                    candidates.push_back(oFace);
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
    if (this->_wasVolumeMesh) {
        auto cmp = _originalMesh.property<Cell, int>("cmp", -1);
        for (auto f : _mesh.faces()) {
            auto oF = _originalMesh.opposite_halfface_handle(_mesh.property(_originalHf, f));
            cmp[_originalMesh.incident_cell(oF)] = _mesh.property(_cmp, f);
        }
    }
}

void StarDecompositionBoundaryLp::apply_current_component() {
    // ASSERT(!_cmpFixV.first.is_valid());

    auto fixVertex = _cmpFixV.first.is_valid() ? _meshVertexMap[_cmpFixV.first] : _mesh.add_vertex_q(_cmpFixV.second);
    for (auto f : _mesh.faces()) {
        if (_mesh.property(_cmp, f) == _cmpIdx) {
            _mesh.delete_face(f, true);
        }
    }

    /*
    for (auto h : _mesh.halfedges()) {
        if (_mesh.is_boundary(h)) {
            auto f = _mesh.add_face(_mesh.from_vertex_handle(h), _mesh.to_vertex_handle(h), fixVertex);
            if (f.is_valid()) {
                _mesh.update_normal_q(f);
                _mesh.property(_cmp, f) = cmpNotSetIdx;
            }
        }
    }
    */

    Mesh& cmpMesh = *_currentCmp;
    fixVertex = _cmpFixV.first.is_valid() ? _cmpFixV.first : cmpMesh.add_vertex_q(_cmpFixV.second);
    for (auto h : cmpMesh.halfedges()) {
        if (cmpMesh.is_boundary(h)) {
            // cmpMesh.add_face(cmpMesh.from_vertex_handle(h), cmpMesh.to_vertex_handle(h), fixVertex);
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
    OpenMesh::VertexHandle borderVertex;
    for (auto fh : _mesh.fh_range(newFace)) {
        auto hv = _mesh.from_vertex_handle(fh);
        auto hvTo = _mesh.to_vertex_handle(fh);
        if (!_cmpVertexMap[hv].is_valid()) {
            auto newVertex = mesh.add_vertex_q(_mesh.data(hv).point_q());
            _cmpVertexMap[hv] = newVertex;
            _meshVertexMap[newVertex] = hv;
            newFaceVertex = hv;
            borderVertex = _cmpVertexMap[hvTo];
        } else if (_currentCmp->find_halfedge(_cmpVertexMap[hvTo], _cmpVertexMap[hv]).is_valid()) {
            existingEdges++;
        } else {
            borderVertex = _cmpVertexMap[hv];
        }

        triangle.push_back(_cmpVertexMap[hv]);
    }

    if (existingEdges == 1 && !newFaceVertex.is_valid()) {
        return false;
    }

    if (!borderVertex.is_valid()) {
        borderVertex = triangle[0];
    }

    auto cmpCenter = has_valid_center(mesh, _currentCmp->data(borderVertex).point_q());
    if (cmpCenter.first == INVALID) {
        if (newFaceVertex.is_valid()) {
            mesh.delete_vertex(_cmpVertexMap[newFaceVertex]);
            _meshVertexMap.erase(_cmpVertexMap[newFaceVertex]);
            _cmpVertexMap.erase(newFaceVertex);
        }
        return false;
    }

    OpenMesh::FaceHandle face = txMesh.add_face(triangle);
    mesh.set_normal_q(face, _mesh.data(newFace).normal_q());
    /*
    for (auto fh : mesh.fh_range(face)) {
        if (mesh.opposite_halfedge_handle(fh).is_boundary()) {
            auto v0 = mesh.to_vertex_handle(fh);
            auto v1 = mesh.from_vertex_handle(fh);
            std::vector<Vector3q> t = { mesh.data(v0).point_q(), mesh.data(v1).point_q(), cmpCenter.second };
            if (triangles_intersect(t, { _meshVertexMap[v0], _meshVertexMap[v1] })) {
                if (newFaceVertex.is_valid()) {
                    _meshVertexMap.erase(_cmpVertexMap[newFaceVertex]);
                    _cmpVertexMap.erase(newFaceVertex);
                    mesh.delete_vertex(newFaceVertex);
                }

                txMesh.revert();
                return false;
            }
        }
    }
    */

    if (cmpCenter.first == VALID_EQUAL) {
        _cmpFixV.first = borderVertex;
    } else {
        _cmpFixV.first = OpenMesh::VertexHandle(-1);
    }
    _cmpFixV.second = cmpCenter.second;
    return true;
}

bool StarDecompositionBoundaryLp::triangles_intersect(std::vector<Vector3q> t, std::vector<OpenMesh::VertexHandle> borderVertices) {
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

std::pair<StarCenterResult, Vector3q> StarDecompositionBoundaryLp::has_valid_center(Mesh& mesh, const Vector3q& point) {
    std::vector<Eigen::Vector3d> positions;
    std::vector<Eigen::Vector3d> normals;
    for (auto face : mesh.faces()) {
        if (face.is_valid() && !face.deleted()) {
            positions.push_back(mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(face))));
            normals.push_back(mesh.normal(face).normalized());
        }
    }

    Eigen::Vector3d p = point.unaryExpr([](mpq_class x) { return x.get_d(); });
    auto result = star_center_close_to(p, positions, normals);
    if (result.first == INVALID) {
        return std::make_pair(result.first, Vector3q::Zero());
    }

    Vector3q newCenter = result.first == VALID_EQUAL ? point : result.second.cast<mpq_class>();
    return std::make_pair(result.first, newCenter);
}

std::pair<bool, Vector3q> StarDecompositionBoundaryLp::has_valid_center(Mesh& mesh) {
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
