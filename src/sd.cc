#include "sd.h"

std::vector<Mesh> sd(Mesh& N) {
    for (auto c : N.cells()) {
        if (N.degenerate_or_inverted(c)) {
            std::cout << "Degenerate or inverted tetrahedron in N" << std::endl;
            return {};
        }
    }

    auto Q = N.property<Vertex, Vector3q>("Q");
    for (auto v : N.vertices()) {
        Q[v] = N.position(v).cast<mpq_class>();
    }


    std::cout << "Decompose" << std::endl;
    // std::vector<Vector3q> centers = decompose(N);
    StarDecompositionBoundary decomposer(N);
    std::vector<Vector3q> centers = decomposer.centers();
    auto cmp = N.property<Cell, int>("cmp");
    int n_cuts = 0;
    for (auto c : N.cells()) {
        int index = cmp[c];
        // TODO: Reenable
        // ASSERT(index >= 0);
        if (n_cuts < index) {
            n_cuts = index;
        }
    }
    std::cout << n_cuts + 1 << " components" << std::endl;
    std::vector<Mesh> components(n_cuts + 1);

    int n_skips = 0;
    int cnt = 0;
    for (int i = 0; i <= n_cuts; i++) {
        Mesh& S = components[i];

        std::vector<Cell> cells;
        for (auto c : N.cells()) {
            if (cmp[c] == i) {
                cells.push_back(c);
                std::vector<Vertex> vertices;
                for (auto cv : N.tet_vertices(c)) {
                    vertices.push_back(S.add_vertex(N.position(cv)));
                }
                S.add_cell(vertices, true);
            }
        }
        N.remove_cells(cells);
    }

    return components;
}
