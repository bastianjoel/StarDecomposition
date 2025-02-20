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

#define FACE_SELECTION_ORDERED 1

const int cmpNotSetIdx = -1;

void StarDecompositionBoundaryLp::run() {
    srand(0);
    std::cout << "Start decomposition" << std::endl;

    _cmpIdx = 0;
    int startOffset = 0;
    while (_mesh.faces_empty() == false) {
        // Add random face
        OpenMesh::FaceHandle h = *(_mesh.faces_begin() + (rand() % _mesh.n_faces()));

        add_component(h);
        if (_mesh.star_center().first) {
            for (auto f : _mesh.faces()) {
                _mesh.property(_cmp, f) = cmpNotSetIdx;
            }
            _cmpMeshes[_cmpIdx] = _mesh;
            break;
        }

        std::vector<OpenMesh::FaceHandle> candidates;
        std::vector<OpenMesh::FaceHandle> candidates2;
        std::set<OpenMesh::FaceHandle> visited;
        for (auto he : _mesh.fh_range(h)) {
            auto oFace = _mesh.opposite_face_handle(he);
            if (oFace.is_valid()) {
                candidates.push_back(oFace);
            }
        }

        while (!candidates.empty() || !candidates2.empty()) {
            if (candidates.empty() && candidates2.empty()) {
                break;
            }

            OpenMesh::FaceHandle nextH;
            if (candidates2.empty()) {
                auto nextPtr = candidates.begin() + (rand() % candidates.size());
                nextH = *nextPtr;
                candidates.erase(nextPtr);
            } else {
                nextH = candidates2.front();
                candidates2.erase(candidates2.begin());
            }

            bool alreadyChecked = _mesh.property(_cmp, nextH) != cmpNotSetIdx || visited.find(nextH) != visited.end();
            if (alreadyChecked || !add_face_to_cmp(*_currentCmp, nextH)) {
                visited.insert(nextH);
                continue;
            }

            std::cout << "." << std::flush;
            visited.clear();
            _viewer->queue_update();
            _mesh.property(_cmp, nextH) = _cmpIdx;
            for (auto he : _mesh.fh_range(nextH)) {
                auto oFace = _mesh.opposite_face_handle(he);
                if (oFace.is_valid() && _mesh.property(_cmp, oFace) == cmpNotSetIdx) {
                    bool between = false;
                    for (auto fh : _mesh.fh_range(oFace)) {
                        auto hv = _mesh.from_vertex_handle(fh);
                        auto hvTo = _mesh.to_vertex_handle(fh);
                        if (fh != _mesh.opposite_halfedge_handle(he) && _cmpVertexMap[hv].is_valid() && _cmpVertexMap[hvTo].is_valid()) {
                            auto h = _currentCmp->find_halfedge(_cmpVertexMap[hv], _cmpVertexMap[hvTo]);
                            if (h.is_valid()) {
                                between = true;
                            }
                            break;
                        }
                    }

                    if (between) {
                        candidates2.push_back(oFace);
                    } else {
                        candidates.push_back(oFace);
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
    }

    std::cout << _cmpIdx + 1 << " components" << std::endl;
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
    for (auto h : cmpMesh.halfedges()) {
        if (cmpMesh.is_boundary(h)) {
            std::vector<OpenMesh::HalfedgeHandle> boundary;
            OpenMesh::HalfedgeHandle next = h;
            do {
                boundary.push_back(next);
                next = cmpMesh.next_halfedge_handle(next);
            } while (next != h);

            std::cout << "Boundary size: " << boundary.size() << std::endl;
            auto cmpFixVertex = cmpMesh.add_vertex_q(_cmpFixV);
            for (auto bH : boundary) {
                cmpMesh.add_face(cmpMesh.from_vertex_handle(bH), cmpMesh.to_vertex_handle(bH), cmpFixVertex);
            }
            // TODO: Move fix vertex
        }
    }

    // TODO: Prevent complex endges
    for (auto h : _mesh.halfedges()) {
        if (_mesh.is_boundary(h)) {
            std::vector<OpenMesh::HalfedgeHandle> boundary;
            OpenMesh::HalfedgeHandle next = h;
            do {
                boundary.push_back(next);
                next = _mesh.next_halfedge_handle(next);
            } while (next != h);

            auto fixVertex = _mesh.add_vertex_q(_cmpFixV);
            for (auto bH : boundary) {
                auto f = _mesh.add_face(_mesh.from_vertex_handle(bH), _mesh.to_vertex_handle(bH), fixVertex);
                if (f.is_valid()) {
                    _mesh.update_normal_q(f);
                    _mesh.property(_cmp, f) = cmpNotSetIdx;
                }

                // TODO: Move fix vertex
            }
        }
    }

    _mesh.delete_isolated_vertices();
    _mesh.garbage_collection();
    _mesh.generate_bvh();
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
    _cmpFixV = _cmpCenter - _mesh.data(startF).normal_q() * 1e-4;

    _mesh.property(_cmp, startF) = _cmpIdx;
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
            OpenMesh::HalfedgeHandle h = _currentCmp->find_halfedge(_cmpVertexMap[hvTo], _cmpVertexMap[hv]);
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
    /*
    if (singleExistingEdge.is_valid() && !newFaceVertex.is_valid()) {
        // Check if one direction consists of a "line" of faces
        auto vExFrom = _meshVertexMap[_currentCmp->from_vertex_handle(singleExistingEdge)];
        auto vExTo = _meshVertexMap[_currentCmp->to_vertex_handle(singleExistingEdge)];
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
                    if (_mesh.property(_cmp, _mesh.opposite_face_handle(fh)) == cmpNotSetIdx) {
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
    }
    */

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
            // auto opposite = _mesh.get_face_in_dir(cmpCenter.second, -meshNormal);

            for (auto h : boundary) {
                auto v0 = mesh.to_vertex_handle(h);
                auto v1 = mesh.from_vertex_handle(h);
                std::vector<Vector3q> t = { mesh.data(v0).point_q(), mesh.data(v1).point_q(), cmpCenter.second - meshNormal * 1e-4 };
                auto intersects = _mesh.triangle_intersects(t, { _meshVertexMap[v0], _meshVertexMap[v1] });
                if (intersects.is_valid()) {
                    cmpCenter.second = cmpCenter.second - meshNormal * 1e-4;
                    valid = false;
                    break;
                }
            }

            if (valid) {
                _cmpCenter = cmpCenter.second;
                _cmpFixV = cmpCenter.second - meshNormal * 1e-4;
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
            _mesh.property(_cmp, f) = _cmpIdx;
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
