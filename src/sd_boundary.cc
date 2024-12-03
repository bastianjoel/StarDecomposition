#include "sd_boundary.h"
#include "lp.h"
#include <map>
#include <queue>
#include <stack>
#include <OpenMesh/Core/IO/MeshIO.hh>

const int cmpNotSetIdx = 0; // TODO: Decrease
StarDecompositionBoundary::StarDecompositionBoundary(Mesh& mesh) : _mesh(mesh), _computed(false) {
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

std::vector<Mesh> StarDecompositionBoundary::components() {
    if (!this->_computed) {
        this->run();
    }

    return std::vector<Mesh>();
}

void StarDecompositionBoundary::run() {
    std::queue<Halfface> queue;
    queue.push(*_surface_hf.begin());
    auto meshCmp = _mesh.property<Cell, int>("cmp", cmpNotSetIdx);

    int cmpIdx = 1; // TODO: Decrease
    while (!queue.empty()) {
        Halfface h = queue.front();
        queue.pop();
        if (_cmp[h] != cmpNotSetIdx) {
            continue;
        }

        std::set<Halfface> cmpHf;
        std::stack<Halfface> candidates;
        candidates.push(h);
        while (!candidates.empty()) {
            auto nextH = candidates.top();
            candidates.pop();
            if (_cmp[nextH] != cmpNotSetIdx) {
                continue;
            }

            if (!add_halfface(cmpHf, nextH)) {
                queue.push(nextH);
                continue;
            }

            meshCmp[_mesh.incident_cell(nextH)] = cmpIdx; // TODO: Remove
            _cmp[nextH] = cmpIdx;
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
        std::cout << "Cmp " << cmpIdx << " amount halffaces: " << cmpHf.size() << std::endl;

        cmpIdx++;
    }

    std::cout << "Num surface halffaces: " << _surface_hf.size() << std::endl;
    std::cout << "Num components: " << cmpIdx << std::endl;
}

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

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
    }

    // Expand boundary to next opposite halfface
    for (auto heh : mesh.halfedges()) {
        if (!mesh.is_boundary(heh)) {
            continue;
        }

        auto v = mesh.to_vertex_handle(heh);
        auto n = mesh.property(normal, mesh.opposite_face_handle(heh));
        // Find cutting boundary face
        for (auto hf : _boundary) {
            auto faceNormal = _normal[hf];
            auto b = faceNormal * _Q[*_mesh.halfface_vertices(hf).first];

            // TODO: Calculate intersection

            // TODO: Check if point is within face bounds
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

            auto newEdge = mesh.find_halfedge(mesh.from_vertex_handle(startHeh), mesh.to_vertex_handle(heh));
            if (newEdge.is_valid() && !newEdge.is_boundary()) {
                continue;
            }

            std::vector<OpenMesh::VertexHandle> triangle;
            triangle.push_back(mesh.from_vertex_handle(startHeh));
            triangle.push_back(mesh.to_vertex_handle(startHeh));
            triangle.push_back(mesh.to_vertex_handle(heh));

            // TODO: Do not add face if it cuts through face in _mesh

            auto face = mesh.add_face(triangle);
            if (!face.is_valid()) {
                throw std::runtime_error("Error: could not add face to component mesh");
                continue;
            }

            Vector3q n = (mesh.property(Q, triangle[1]) - mesh.property(Q, triangle[0])).cross(mesh.property(Q, triangle[2]) - mesh.property(Q, triangle[0]));
            mesh.property(normal, face) = n;
            gotToEnd = false;
            addedFaces++;
            break;
        }
    }

    std::cout << "Added faces: " << addedFaces << std::endl;

    std::vector<Eigen::Vector3d> positions;
    std::vector<Eigen::Vector3d> normals;
    for (auto face : mesh.faces()) {
        positions.push_back(mesh.property(Q, mesh.to_vertex_handle(mesh.halfedge_handle(face))).unaryExpr([](mpq_class x) { return x.get_d(); }));
        normals.push_back(mesh.property(normal, face).unaryExpr([](mpq_class x) { return x.get_d(); }).normalized());
    }

    Vector3q new_center = kernel_chebyshev_center(positions, normals).cast<mpq_class>();
    bool new_center_valid = true;
    for (auto face : mesh.faces()) {
        std::vector<Vector3q> triangle;
        for (auto fv : mesh.fv_range(face)) {
            triangle.push_back(mesh.property(Q, fv));
        }
        Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
        if ((triangle[0] - new_center).dot(n) <= 0) {
            new_center_valid = false;
            break;
        }
    }

    if (!new_center_valid) {
        if (addedFaces > 10 && !OpenMesh::IO::write_mesh(mesh, "invalid_center.off")) {
            std::cerr << "write error\n";
            exit(1);
        }

        cmpHf.erase(hf);
    } else {
        if (addedFaces > 10 && !OpenMesh::IO::write_mesh(mesh, "valid_center.off")) {
            std::cerr << "write error\n";
            exit(1);
        }
    }

    return new_center_valid;
}
