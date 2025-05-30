#include "sd_boundary.h"
#include "retet.h"
#include "sd.h"

StarDecompositionBoundary::StarDecompositionBoundary(Mesh& m, int seed) : _mesh(m), seed(seed) {
    _mesh.add_property(_selected);
    _mesh.add_property(_origBound);
    for (auto v : _mesh.vertices()) {
        _mesh.data(v).set_point_q(_mesh.point(v).cast<mpq_class>());
    }

    for (auto f : _mesh.faces()) {
        _mesh.update_normal_q(f);
        _mesh.property(_selected, f) = false;
        _mesh.property(_origBound, f) = true;
    }

    _mesh.generate_bvh();

#ifdef SAVE_DEBUG_MESHES
    if (!OpenMesh::IO::write_mesh(_mesh, "debug/original.obj")) {
        std::cerr << "write error\n";
        exit(1);
    }
#endif
}

StarDecompositionBoundary::StarDecompositionBoundary(VolumeMesh& mesh, int seed) : seed(seed) {
    _mesh.add_property(_selected);
    _mesh.add_property(_origBound);

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
        _mesh.property(_selected, face) = false;
        _mesh.property(_origBound, face) = true;
    }

    _mesh.delete_isolated_vertices();
    _mesh.garbage_collection();
    _mesh.generate_bvh();

#ifdef SAVE_DEBUG_MESHES
    if (!OpenMesh::IO::write_mesh(_mesh, "debug/original.obj")) {
        std::cerr << "write error\n";
        exit(1);
    }
#endif
}

std::vector<Vector3q> StarDecompositionBoundary::centers() {
    if (!this->_computed) {
        this->start();
    }

    return std::vector<Vector3q>();
}

std::vector<VolumeMesh> StarDecompositionBoundary::components() {
    if (!this->_computed) {
        this->start();
    }

    return _components;
}

void StarDecompositionBoundary::add_component(const Mesh& mesh) {
#ifdef SAVE_DEBUG_MESHES
    char buffer[100];
    std::sprintf(buffer, "debug/cmp_%zu.obj", _components.size());
    if (!OpenMesh::IO::write_mesh(mesh, buffer)) {
        std::cerr << "write error\n";
        exit(1);
    }

    std::sprintf(buffer, "debug/original_%zu.obj", _components.size());
    if (!OpenMesh::IO::write_mesh(_mesh, buffer)) {
        std::cerr << "write error\n";
        exit(1);
    }
#endif

    VolumeMesh m;
    if (!retetrahedrize(mesh, m)) {
        std::cerr << "Could not tetrahedrize component" << std::endl;
        exit(1);
    } else {
        _components.push_back(m);
    }
}

bool StarDecompositionBoundary::triangle_intersects(const OpenMesh::HalfedgeHandle& h, const Vector3q& p) {
    auto v0 = _nextComponent.to_vertex_handle(h);
    auto v1 = _nextComponent.from_vertex_handle(h);

    std::vector<Vector3q> t = { _nextComponent.data(v0).point_q(), _nextComponent.data(v1).point_q(), p };
    return _mesh.triangle_intersects(t, { _meshVertexMap[v0], _meshVertexMap[v1] }).is_valid();
}

bool StarDecompositionBoundary::is_valid_with(const MeshBoundary& b, const OpenMesh::FaceHandle& f, const Vector3q& c) {
    if ((_mesh.data(*_mesh.fv_begin(f)).point_q() - c).dot(_mesh.data(f).normal_q()) <= 0) {
        return false;
    }

    for (auto fh : _mesh.fh_range(f)) {
        auto hv = _mesh.from_vertex_handle(fh);
        auto hvTo = _mesh.to_vertex_handle(fh);
        if (!_nextComponent.find_halfedge(_cmpVertexMap[hv], _cmpVertexMap[hvTo]).is_valid()) {
            auto v0 = _mesh.data(hvTo).point_q();
            auto v1 = _mesh.data(hv).point_q();
            auto n = (v1 - v0).cross(b.get_fix_vertex() - v0);
            if ((v0 - c).dot(n) <= 0) {
                return false;
            }

            std::vector<Vector3q> t = { v0, v1, b.get_fix_vertex() };
            if (_mesh.triangle_intersects(t, { hv, hvTo }).is_valid()) {
                return false;
            }
        }
    }

    return true;
}

std::optional<Vector3q> StarDecompositionBoundary::is_next_component_valid() {
    std::vector<Eigen::Vector3d> positions;
    std::vector<Eigen::Vector3d> normals;
    for (auto face : _nextComponent.faces()) {
        positions.push_back(_nextComponent.point(_nextComponent.to_vertex_handle(_nextComponent.halfedge_handle(face))));
        normals.push_back(_nextComponent.data(face).normal_normalized());
    }

    for (auto b : _boundaries) {
        for (auto n : b.get_normals()) {
            positions.push_back(b.get_fix_vertex().unaryExpr([](mpq_class x) { return x.get_d(); }));
            normals.push_back(n);
        }
    }

    auto newCenter = kernel_chebyshev_center(positions, normals);
    for (int i = 0; i < positions.size(); i++) {
        if ((positions[i] - newCenter).dot(normals[i]) <= 0) {
            return std::nullopt;
        }
    }

    return newCenter.cast<mpq_class>();
}

void StarDecompositionBoundary::fallback(const Mesh& mesh) {
    auto components = sd(_mesh, "tet", seed);
    for (auto m : components) {
        _components.push_back(m);
    }
}

void StarDecompositionBoundary::start() {
    int startOffset = 0;
    int resets = 0;
    int bestSinceReset = 0;
    while (_mesh.faces_empty() == false) {
        // Check if mesh is already star shaped
        if (_mesh.star_center().first) {
            add_component(_mesh);
            break;
        }

        // Add random face
        OpenMesh::FaceHandle h = *(_mesh.faces_begin() + (rand() % _mesh.n_faces()));
        while (!_mesh.property(_origBound, h) && startOffset < _mesh.n_faces()) {
            h = *(_mesh.faces_begin() + (startOffset));
            startOffset++;
        }

        if (startOffset == _mesh.n_faces()) {
            fallback(_mesh);
            break;
        }

        _cmpVertexMap = std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle>();
        _meshVertexMap = std::map<OpenMesh::VertexHandle, OpenMesh::VertexHandle>();
        _boundaries = std::vector<MeshBoundary>();
        _nextComponent = Mesh();
        _nextComponentFaces = 0;
        _boundaries.push_back(MeshBoundary(_nextComponent));
        init_component(_nextComponent, h);
        for (auto f :_nextComponent.faces()) {
            _nextComponentFaces++;
        }

        std::vector<OpenMesh::FaceHandle> candidates;
        std::vector<OpenMesh::FaceHandle> candidates2;
        std::vector<OpenMesh::FaceHandle> failedRecoverable;
        std::set<OpenMesh::FaceHandle> visited;
        for (auto he : _mesh.fh_range(h)) {
            auto oFace = _mesh.opposite_face_handle(he);
            if (oFace.is_valid()) {
                candidates.push_back(oFace);
            }
        }

        while (!candidates.empty() || !candidates2.empty()) {
            OpenMesh::FaceHandle nextH;
            if (candidates2.empty()) {
                auto nextPtr = candidates.begin() + (rand() % candidates.size());
                nextH = *nextPtr;
                candidates.erase(nextPtr);
            } else {
                auto nextPtr = candidates2.begin() + (rand() % candidates2.size());
                nextH = *nextPtr;
                candidates2.erase(nextPtr);
            }

            bool alreadyChecked = _mesh.property(_selected, nextH) || visited.find(nextH) != visited.end();
            if (alreadyChecked) {
                visited.insert(nextH);
                continue;
            }

            int addStat = add_face_to_cmp(_nextComponent, nextH);
            if (addStat != 0) {
                visited.insert(nextH);
                if (addStat <= 2) {
                    failedRecoverable.push_back(nextH);
                }
                continue;
            }

            if (_recheckFailed) {
                for (auto f : failedRecoverable) {
                    if (_mesh.property(_selected, f)) {
                        continue;
                    }

                    candidates.push_back(f);
                }
                failedRecoverable.clear();
                _recheckFailed = false;
            }
#ifdef GUI
            _viewer->queue_update();
#endif
            visited.clear();
            _nextComponentFaces++;
            _mesh.property(_selected, nextH) = true;
            for (auto he : _mesh.fh_range(nextH)) {
                auto oFace = _mesh.opposite_face_handle(he);
                if (oFace.is_valid() && !_mesh.property(_selected, oFace)) {
                    bool between = false;
                    for (auto fh : _mesh.fh_range(oFace)) {
                        auto hv = _mesh.from_vertex_handle(fh);
                        auto hvTo = _mesh.to_vertex_handle(fh);
                        if (fh != _mesh.opposite_halfedge_handle(he) && _cmpVertexMap[hv].is_valid() && _cmpVertexMap[hvTo].is_valid()) {
                            auto h = _nextComponent.find_halfedge(_cmpVertexMap[hv], _cmpVertexMap[hvTo]);
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

        int numCmpFacesTotal = _nextComponent.faces().to_vector().size();
        int numCmpFaces = 0;
        for (auto f : _mesh.faces()) {
            if (_mesh.property(_selected, f)) {
                numCmpFaces++;
            }
        }
        int amount = (-numCmpFaces + (numCmpFacesTotal - numCmpFaces));
        // if (amount < -10 || (resets > 3 && amount <= 0.7 * bestSinceReset && amount < 0)) {
        if (amount < -50 || (resets > 8 && amount <= 0.7 * bestSinceReset && amount < -10)) {
            std::cout << "Cmp " << _components.size() << " amount faces: " << _nextComponent.faces().to_vector().size() << " full mesh size change: " << (-numCmpFaces + (numCmpFacesTotal - numCmpFaces));
            std::cout << " resets: " << resets << " full mesh size: " << _mesh.n_faces() << std::endl;
            finalize_component(_nextComponent);
            add_component(_nextComponent);
            startOffset = 0;
            resets = 0;
            bestSinceReset = 0;
#ifdef GUI
            _viewer->clear_extras();
#endif
        } else if (resets > 10) {
            fallback(_mesh);
            break;
        } else {
            for (auto f : _mesh.faces()) {
                if (_mesh.property(_selected, f)) {
                    _mesh.property(_selected, f) = false;
                }
            }
        
            startOffset++;
            resets++;

            if (amount < bestSinceReset) {
                bestSinceReset = (-numCmpFaces + (numCmpFacesTotal - numCmpFaces));
            }
        }
    }

    this->_computed = true;
}
