#include "OpenMesh/Core/Utils/Property.hh"
#include "sd.h"
#include <stack>

// Creates a star decomposition of the mesh.
// Uses the outer faces to create components by adding neighboring faces
// that are visible from a point within the component.
std::vector<Vector3q> decompose_by_corners(Mesh& mesh) {
    std::set<Halfface> surface = get_surface_halffaces(mesh);
    std::queue<Halfface> queue;
    queue.push(*surface.begin());
    const int cmpNotSetIdx = 0; // TODO: Decrease
    auto meshCmp = mesh.property<Cell, int>("cmp", cmpNotSetIdx);
    auto meshHfCmp = mesh.property<Halfface, int>("cmp", cmpNotSetIdx);

    auto Q = mesh.property<Vertex, Vector3q>("Q");
    auto normal = mesh.property<Halfface, Vector3q>();
    for (auto hf : surface) {
        Halfface h = mesh.halfface_handle(hf.face_handle(), 0);
        std::vector<Vector3q> triangle;
        for (auto hv : mesh.halfface_vertices(h)) {
            triangle.push_back(Q[hv]);
        }
        Vector3q n = (triangle[1] - triangle[0]).cross(triangle[2] - triangle[0]);
        normal[h] = n;
        normal[mesh.opposite_halfface_handle(h)] = -n;
    }

    int cmpIdx = 1; // TODO: Decrease
    while (!queue.empty()) {
        Halfface h = queue.front();
        queue.pop();
        if (meshHfCmp[h] != cmpNotSetIdx) {
            continue;
        }

        Mesh cmpMesh;
        auto cmpVertexMap = mesh.property<Vertex, Vertex>();
        std::stack<Halfface> candidates;
        candidates.push(h);
        while (!candidates.empty()) {
            auto nextH = candidates.top();
            candidates.pop();
            if (meshHfCmp[nextH] != cmpNotSetIdx) {
                continue;
            }

            if (!can_add_halfface(mesh, cmpMesh, nextH, cmpVertexMap)) {
                queue.push(nextH);
                continue;
            }

            meshCmp[mesh.incident_cell(nextH)] = cmpIdx; // TODO: Remove
            meshHfCmp[nextH] = cmpIdx;
            for (auto hv : mesh.halfface_vertices(nextH)) {
                for (auto vh : mesh.vertex_halffaces(hv)) {
                    // Skip if cell is not a surface cell or starting halfface
                    if (surface.find(vh) == surface.end() || nextH == vh) {
                        continue;
                    }

                    if (meshHfCmp[vh] != cmpNotSetIdx) {
                        continue;
                    }

                    candidates.push(vh);
                }
            }
        }
        std::cout << "Cmp " << cmpIdx << " amount halffaces: " << cmpMesh.n_halffaces() << std::endl;

        cmpIdx++;
    }

    std::cout << "Num surface halffaces: " << surface.size() << std::endl;
    std::cout << "Num components: " << cmpIdx << std::endl;

    return std::vector<Vector3q>();
}

std::set<Halfface> get_surface_halffaces(Mesh& mesh) {
    std::set<Halfface> surface;
    for (auto h : mesh.halffaces()) {
        if (mesh.is_boundary(h)) {
            surface.insert(h.opposite_handle());
        }
    }

    return surface;
}

bool can_add_halfface(Mesh& parentMesh, Mesh& cmpMesh, Halfface& hf, Property<Vertex, Vertex> vertexMap) {
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
