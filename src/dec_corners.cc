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

    auto Q = mesh.property<Vertex, Vector3q>("Q");
    auto normal = mesh.property<Halfface, Vector3q>();
    for (auto f : mesh.faces()) {
        Halfface h = mesh.halfface_handle(f, 0);
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
        Cell hc = mesh.incident_cell(h);
        queue.pop();
        if (meshCmp[hc] != cmpNotSetIdx) {
            continue;
        }

        meshCmp[hc] = cmpIdx;
        std::stack<Halfface> candidates;
        candidates.push(h);
        while (!candidates.empty()) {
            auto nextH = candidates.top();
            candidates.pop();
            // TODO: Check if `nextH` can and should be added
            int amount = (normal[nextH][0] == normal[h][0]) + (normal[nextH][1] == normal[h][1]) + (normal[nextH][2] == normal[h][2]);
            if (amount < 2) {
                queue.push(nextH);
                continue;
            }

            meshCmp[mesh.incident_cell(nextH)] = cmpIdx;
            for (auto hv : mesh.halfface_vertices(nextH)) {
                for (auto vh : mesh.vertex_halffaces(hv)) {
                    // Skip if cell is not a surface cell or starting halfface
                    if (surface.find(vh) == surface.end() || nextH == vh) {
                        continue;
                    }

                    if (meshCmp[mesh.incident_cell(vh)] != cmpNotSetIdx) {
                        continue;
                    }

                    candidates.push(vh);
                }
            }
        }

        cmpIdx++;
    }

    std::cout << "Num surface halffaces: " << surface.size() << std::endl;
    std::cout << "Num components: " << cmpIdx << std::endl;

    return std::vector<Vector3q>();
}

std::set<Halfface> get_surface_halffaces(Mesh& mesh) {
    std::set<Halfface> surface;
    for (auto f : mesh.faces()) {
        for (auto h : mesh.face_halffaces(f)) {
            if (mesh.is_boundary(h)) {
                surface.insert(h.opposite_handle());
            }
        }
    }

    return surface;
}
