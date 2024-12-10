#include "sd_boundary.h"
#include "lp.h"
#include <cstdio>
#include <map>
#include <queue>
#include <stack>
#include <OpenMesh/Core/IO/MeshIO.hh>

const int cmpNotSetIdx = -1;
StarDecompositionBoundary::StarDecompositionBoundary(VolumeMesh& mesh) : _mesh(mesh), _computed(false) {
    _cmp = mesh.property<Halfface, int>("cmp", cmpNotSetIdx);

    _Q = mesh.property<Vertex, Vector3q>("Q");
    _normal = mesh.property<Halfface, Vector3q>();
    _boundary = _mesh.boundary_halffaces();
    for (auto hf : _boundary) {
        Halfface h = mesh.halfface_handle(hf.face_handle(), 0);
        std::vector<Vector3q> triangle;
        for (auto hv : mesh.halfface_vertices(h)) {
            triangle.push_back(_Q[hv]);
        }
        Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
        _normal[h] = n;
        _normal[mesh.opposite_halfface_handle(h)] = -n;
        _surface_hf.insert(h);
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
    std::queue<Halfface> queue;
    queue.push(*_surface_hf.begin());
    auto meshCmp = _mesh.property<Cell, int>("cmp", cmpNotSetIdx);

    _cmpIdx = 0;
    while (!queue.empty()) {
        Halfface h = queue.front();
        queue.pop();
        if (_cmp[h] != cmpNotSetIdx) {
            continue;
        }

        add_component(h);

        std::set<Halfface> cmpHf;
        std::stack<Halfface> candidates;
        for (auto he : _mesh.halfface_edges(h)) {
            for (auto eh : _mesh.edge_halffaces(he)) {
                if (_surface_hf.find(eh) != _surface_hf.end() && _cmp[eh] == cmpNotSetIdx) {
                    candidates.push(eh);
                }
            }
        }

        while (!candidates.empty()) {
            auto nextH = candidates.top();
            candidates.pop();
            if (_cmp[nextH] != cmpNotSetIdx) {
                continue;
            }

            if (!add_hf_to_cmp(_cmpIdx, nextH)) {
                queue.push(nextH);
                continue;
            }

            meshCmp[_mesh.incident_cell(nextH)] = _cmpIdx; // TODO: Remove
            _cmp[nextH] = _cmpIdx;
            for (auto he : _mesh.halfface_edges(nextH)) {
                for (auto eh : _mesh.edge_halffaces(he)) {
                    // Skip if cell is not a surface cell or starting halfface
                    if (_surface_hf.find(eh) == _surface_hf.end() || nextH == eh) {
                        continue;
                    }

                    if (_cmp[eh] != cmpNotSetIdx) {
                        continue;
                    }

                    candidates.push(eh);
                }
            }
        }
        std::cout << "Cmp " << _cmpIdx << " amount faces: " << cmpHf.size() << std::endl;
    }

    std::cout << "Num surface halffaces: " << _surface_hf.size() << std::endl;
    std::cout << "Num components: " << _cmpIdx << std::endl;
}

MyMesh StarDecompositionBoundary::add_component(Halfface& startHf) {
    _cmpVertexMap = std::map<Vertex, OpenMesh::SmartVertexHandle>();

    MyMesh mesh;
    std::vector<OpenMesh::VertexHandle> newHfVertices;
    for (auto hv : _mesh.halfface_vertices(startHf)) {
        auto newVertex = mesh.add_vertex(MyMesh::Point(_Q[hv][0].get_d(), _Q[hv][1].get_d(), _Q[hv][2].get_d()));
        mesh.data(newVertex).set_q(_Q[hv]);
        _cmpVertexMap[hv] = newVertex;
        newHfVertices.push_back(newVertex);
    }
    auto face = mesh.add_face(newHfVertices);
    mesh.set_normal(face, _normal[startHf].unaryExpr([](mpq_class x) { return x.get_d(); }));

    _cmpMeshes.push_back(mesh);
    _cmpIdx = _cmpMeshes.size() - 1;

    // TODO: Remove setting cell component property
    auto meshCmp = _mesh.property<Cell, int>("cmp", cmpNotSetIdx);
    meshCmp[_mesh.incident_cell(startHf)] = _cmpIdx;

    _cmp[startHf] = _cmpIdx;
    return mesh;
}

bool StarDecompositionBoundary::add_hf_to_cmp(int cmp, Halfface& hf) {
    std::vector<Vertex> newHfVertices;
    std::vector<OpenMesh::SmartVertexHandle> triangle;
    for (auto hv : _mesh.halfface_vertices(hf)) {
        if (!_cmpVertexMap[hv].is_valid()) {
            auto newVertex = _cmpMeshes[cmp].add_vertex(MyMesh::Point(_Q[hv][0].get_d(), _Q[hv][1].get_d(), _Q[hv][2].get_d()));
            _cmpMeshes[cmp].data(newVertex).set_q(_Q[hv]);
            _cmpVertexMap[hv] = newVertex;
            newHfVertices.push_back(hv);
        }

        triangle.push_back(_cmpVertexMap[hv]);
    }

    auto face = _cmpMeshes[cmp].add_face(triangle);

    bool hasValidCenter = has_valid_center(_cmpMeshes[cmp]);
    if (!hasValidCenter) {
        _cmpMeshes[cmp].delete_face(face, true);
        for (auto v : newHfVertices) {
            _cmpVertexMap.erase(v);
        }
    }

    return hasValidCenter;
}

bool StarDecompositionBoundary::add_halfface(std::set<Halfface>& cmpHf, Halfface& hf) {
    cmpHf.insert(hf);

    if (cmpHf.size() == 1) {
        return true;
    }

    MyMesh mesh;
    OpenMesh::FPropHandleT<Vector3q> normal;
    OpenMesh::VPropHandleT<Vector3q> Q;
    mesh.add_property(normal);
    mesh.add_property(Q);
    std::map<Vertex, OpenMesh::SmartVertexHandle> vertexMap;
    for (auto h : cmpHf) {
        std::vector<OpenMesh::SmartVertexHandle> newHfVertices;
        for (auto hv : _mesh.halfface_vertices(h)) {
            if (!vertexMap[hv].is_valid()) {
                auto newVertex = mesh.add_vertex(MyMesh::Point(_Q[hv][0].get_d(), _Q[hv][1].get_d(), _Q[hv][2].get_d()));
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
            auto newV = mesh.add_vertex(MyMesh::Point(newVPos[0].get_d(), newVPos[1].get_d(), newVPos[2].get_d()));
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

bool StarDecompositionBoundary::has_valid_center(MyMesh mesh) {
    std::vector<Eigen::Vector3d> positions;
    std::vector<Eigen::Vector3d> normals;
    for (auto face : mesh.faces()) {
        positions.push_back(mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(face))));
        normals.push_back(mesh.normal(face).normalized());
    }

    Vector3q new_center = kernel_chebyshev_center(positions, normals).cast<mpq_class>();
    for (auto face : mesh.faces()) {
        std::vector<Vector3q> triangle;
        for (auto fv : mesh.fv_range(face)) {
            triangle.push_back(mesh.data(fv).Q());
        }
        Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
        if ((triangle[0] - new_center).dot(n) <= 0) {
            return false;
        }
    }

    return true;
}
