#include "sd.h"
#include "retet.h"
#include "sd_boundary_lp.h"

#ifdef GUI
#include "viewer.h"

VolumeMesh viewer_mesh;
Viewer viewer(viewer_mesh);

bool show_coordinate_system = false;

void callback() {
    if (ImGui::IsKeyPressed(ImGuiKey_B)) {
        viewer.show_boundary_ = !viewer.show_boundary_;
        viewer.update();
    }
    if (ImGui::IsKeyPressed(ImGuiKey_M)) {
        viewer.show_colored_ = !viewer.show_colored_;
        viewer.update();
    }
    if (ImGui::IsKeyPressed(ImGuiKey_Q)) {
        viewer.reset();
    }
    if (ImGui::IsKeyPressed(ImGuiKey_X)) {
        show_coordinate_system = !show_coordinate_system;
        viewer.clear_extras();
        if (show_coordinate_system) {
            viewer.add_coordinate_system();
        }
        viewer.update();
    }
}
#endif

/*
#ifdef GUI
        auto clr = viewer_mesh.property<Cell, Eigen::Vector3d>("clr", {1, 1, 1});
        VolumeMesh& Mi = N;
        Eigen::Vector3d color = {1, 1, 1};
        auto vmap = Mi.property<Vertex, Vertex>();
        for (auto c : Mi.cells()) {
            std::vector<Vertex> vertices;
            for (auto cv : Mi.tet_vertices(c)) {
                if (!viewer_mesh.is_valid(vmap[cv])) {
                    Vertex v = viewer_mesh.add_vertex(Mi.position(cv));
                    vmap[cv] = v;
                }
                vertices.push_back(vmap[cv]);
            }
            clr[viewer_mesh.add_cell(vertices, true)] = color;
        }
        components = sd(N, algorithm, &viewer);
        viewer_thread.join();
#endif
*/

std::vector<VolumeMesh> sd(VolumeMesh& N, std::string algorithm, int seed) {
    for (auto c : N.cells()) {
        if (N.degenerate_or_inverted(c)) {
            std::cout << "Degenerate or inverted tetrahedron in N" << std::endl;
            return {};
        }
    }

    bool qNeedInit = !N.property_exists<Vector3q, typename Vertex::EntityTag>("Q");
    auto Q = N.property<Vertex, Vector3q>("Q");
    if (qNeedInit) {
        for (auto v : N.vertices()) {
            Q[v] = N.position(v).cast<mpq_class>();
        }
    }

#ifdef GUI
    auto clr = viewer_mesh.property<Cell, Eigen::Vector3d>("clr", {1, 1, 1});
    VolumeMesh& Mi = N;
    Eigen::Vector3d color = {1, 1, 1};
    auto vmap = Mi.property<Vertex, Vertex>();
    for (auto c : Mi.cells()) {
        std::vector<Vertex> vertices;
        for (auto cv : Mi.tet_vertices(c)) {
            if (!viewer_mesh.is_valid(vmap[cv])) {
                Vertex v = viewer_mesh.add_vertex(Mi.position(cv));
                vmap[cv] = v;
            }
            vertices.push_back(vmap[cv]);
        }
        clr[viewer_mesh.add_cell(vertices, true)] = color;
    }
#endif

    std::vector<Vector3q> centers;
    std::cout << "Decompose using ";
    if (algorithm == "tet") {
        std::cout << "tetrahedron" << std::endl;
        centers = decompose(N, seed);
    } else if (algorithm == "boundary") {
        std::cout << "boundary" << std::endl;
        StarDecompositionBoundaryChebyshev decomposer(N, seed);
#ifdef GUI
        decomposer.set_viewer(&viewer);
        std::thread viewer_thread([]() {
            viewer.start(callback, "..", "Star Decomposition Maps");
        });
#endif
        centers = decomposer.centers();
#ifdef GUI
        viewer_thread.join();
#endif
        return decomposer.components();
    } else {
        std::cout << "boundary (lp center move)" << std::endl;
        StarDecompositionBoundaryLp decomposer(N, seed);
#ifdef GUI
        decomposer.set_viewer(&viewer);
        std::thread viewer_thread([]() {
            viewer.start(callback, "..", "Star Decomposition Maps");
        });
#endif
        centers = decomposer.centers();
#ifdef GUI
        viewer_thread.join();
#endif
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
    std::vector<VolumeMesh> components(n_cuts + 1);

    int n_skips = 0;
    int cnt = 0;
    for (int i = 0; i <= n_cuts; i++) {
        VolumeMesh& S = components[i];

        std::map<Vertex, Vertex> vmap;
        std::vector<Cell> cells;
        for (auto c : N.cells()) {
            if (cmp[c] == i) {
                cells.push_back(c);
                std::vector<Vertex> vertices;
                for (auto cv : N.tet_vertices(c)) {
                    if (vmap.find(cv) == vmap.end()) {
                        vmap[cv] = S.add_vertex(N.position(cv));
                    }

                    vertices.push_back(vmap[cv]);
                }
                S.add_cell(vertices, true);
            }
        }
        N.remove_cells(cells);
    }

    return components;
}

std::vector<VolumeMesh> sd(Mesh& N, std::string algorithm, int seed) {
#ifdef GUI
    VolumeMesh M;
    if (!retetrahedrize(N, M)) {
        return {};
    }

    return sd(M, algorithm, seed);
#endif
#ifndef GUI
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
        StarDecompositionBoundaryChebyshev decomposer(N, seed);
        std::vector<Vector3q> centers = decomposer.centers();
        return decomposer.components();
    }

    std::cout << "boundary (lp center move)" << std::endl;
    StarDecompositionBoundaryLp decomposer(N, seed);
    std::vector<Vector3q> centers = decomposer.centers();

    return decomposer.components();
#endif
}
