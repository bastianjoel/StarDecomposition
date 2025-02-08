#include "sd.h"
#include "retet.h"
#include "sd_boundary_lp.h"

std::vector<VolumeMesh> sd(VolumeMesh& N, std::string algorithm) {
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

    std::vector<Vector3q> centers;
    std::cout << "Decompose using ";
    if (algorithm == "tet") {
        std::cout << "tetrahedron" << std::endl;
        centers = decompose(N);
    } else if (algorithm == "boundary") {
        std::cout << "boundary" << std::endl;
        StarDecompositionBoundary decomposer(N);
        centers = decomposer.centers();
        return decomposer.components();
    } else {
        std::cout << "boundary (lp center move)" << std::endl;
        StarDecompositionBoundaryLp decomposer(N);
        centers = decomposer.centers();
        return decomposer.components();
    }

    auto cmp = N.property<Cell, int>("cmp");
    int n_cuts = 0;
    for (auto c : N.cells()) {
        int index = cmp[c];
        // ASSERT(index >= 0);
        if (n_cuts < index) {
            n_cuts = index;
        }
    }
    std::cout << n_cuts + 1 << " components" << std::endl;
    std::vector<VolumeMesh> components(n_cuts + 1);

    int n_skips = 0;
    int cnt = 0;
    for (int i = 0; i <= n_cuts; i++) {
        VolumeMesh& S = components[i];

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

std::vector<VolumeMesh> sd(Mesh& N, std::string algorithm) {
    if (algorithm == "tet") {
        VolumeMesh M;
        if (!retetrahedrize(N, M)) {
            return {};
        }

        return sd(M, algorithm);
    }

    std::cout << "Decompose using ";
    if (algorithm == "boundary") {
        std::cout << "boundary" << std::endl;
        StarDecompositionBoundary decomposer(N);
        std::vector<Vector3q> centers = decomposer.centers();
        return decomposer.components();
    }

    std::cout << "boundary (lp center move)" << std::endl;
    StarDecompositionBoundaryLp decomposer(N);
    std::vector<Vector3q> centers = decomposer.centers();

    return decomposer.components();
}
