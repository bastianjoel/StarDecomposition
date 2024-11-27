#include "sd_boundary.h"
#include <queue>
#include <stack>

const int cmpNotSetIdx = 0; // TODO: Decrease
StarDecompositionBoundary::StarDecompositionBoundary(Mesh& mesh) : _mesh(mesh), _computed(false) {
    _cmp = mesh.property<Halfface, int>("cmp", cmpNotSetIdx);

    auto Q = mesh.property<Vertex, Vector3q>("Q");
    _normal = mesh.property<Halfface, Vector3q>();
    for (auto hf : _mesh.boundary_halffaces()) {
        Halfface h = mesh.halfface_handle(hf.face_handle(), 0);
        std::vector<Vector3q> triangle;
        for (auto hv : mesh.halfface_vertices(h)) {
            triangle.push_back(Q[hv]);
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

        Mesh cmpMesh;
        auto cmpVertexMap = _mesh.property<Vertex, Vertex>();
        std::stack<Halfface> candidates;
        candidates.push(h);
        while (!candidates.empty()) {
            auto nextH = candidates.top();
            candidates.pop();
            if (_cmp[nextH] != cmpNotSetIdx) {
                continue;
            }

            if (!can_add_halfface(cmpMesh, nextH, cmpVertexMap)) {
                queue.push(nextH);
                continue;
            }

            meshCmp[_mesh.incident_cell(nextH)] = cmpIdx; // TODO: Remove
            _cmp[nextH] = cmpIdx;
            for (auto hv : _mesh.halfface_vertices(nextH)) {
                for (auto vh : _mesh.vertex_halffaces(hv)) {
                    // Skip if cell is not a surface cell or starting halfface
                    if (_surface_hf.find(vh) == _surface_hf.end() || nextH == vh) {
                        continue;
                    }

                    if (_cmp[vh] != cmpNotSetIdx) {
                        continue;
                    }

                    candidates.push(vh);
                }
            }
        }
        std::cout << "Cmp " << cmpIdx << " amount halffaces: " << cmpMesh.n_halffaces() << std::endl;

        cmpIdx++;
    }

    std::cout << "Num surface halffaces: " << _surface_hf.size() << std::endl;
    std::cout << "Num components: " << cmpIdx << std::endl;
}

bool StarDecompositionBoundary::can_add_halfface(Mesh& cmpMesh, Halfface& hf, Property<Vertex, Vertex> vertexMap) {
    /*
    std::vector<Vertex> vertices;
    std::vector<Vertex> verticesAdded;
    std::vector<Vertex> connecting;
    for (auto hv : parentMesh.halfface_vertices(hf)) {
        if (!vertexMap[hv].is_valid()) {
            auto newVertex = cmpMesh.add_vertex(parentMesh.position(hv));
            vertexMap[hv] = newVertex;
            vertices.push_back(newVertex);
            verticesAdded.push_back(newVertex);
        } else {
            vertices.push_back(vertexMap[hv]);
            connecting.push_back(vertexMap[hv]);
        }
    }
    auto newHf = cmpMesh.add_halfface(vertices[0], vertices[1], vertices[2], true);
    */

    /* 
        int amount = (normal[nextH][0] == normal[h][0]) + (normal[nextH][1] == normal[h][1]) + (normal[nextH][2] == normal[h][2]);
        if (amount < 2) {
            queue.push(nextH);
            continue;
        }
    */

    return true;
}
