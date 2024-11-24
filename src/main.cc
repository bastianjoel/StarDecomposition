#include "retet.h"
#include "sd.h"

#ifdef GUI
#include "viewer.h"

Mesh viewer_mesh;
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

void print_usage(std::string arg0, std::string path) {
    std::cout << "Usage: " << arg0 << " <N> [-o <out_dir>]" << std::endl;
    std::cout << "  e.g. " << arg0 << " " << path << SEP << "meshes" << SEP << "thumb.vtk" << std::endl;
}

int main(int argc, char** argv) {
    std::string arg0 = std::string(argv[0]);
    std::string path = arg0.substr(0, arg0.rfind(SEP) + 1) + "..";
    std::string N_filename = "";
    std::string out_dir = "";
    for (int i = 1; i < argc; i++) {
        std::string arg(argv[i]);
        if (arg.at(0) == '-') {
            for (int j = 1; j < arg.size(); j++) {
                char c = arg.at(j);
                switch (c) {
                case 'h':
                    print_usage(arg0, path);
                    return 0;
                case 'o':
                    i++;
                    if (argc <= i) {
                        std::cout << "Output directory missing" << std::endl;
                        print_usage(arg0, path);
                        return 0;
                    }
                    out_dir = std::string(argv[i]);
                    break;
                default:
                    std::cout << "Invalid option \"" << c << "\"" << std::endl;
                    print_usage(arg0, path);
                    return 0;
                }
            }
        } else if (N_filename == "") {
            N_filename = arg;
        }
    }
    if (N_filename == "") {
        N_filename = path + SEP + "meshes" + SEP + "thumb.vtk";
    }
    int orientation = 0;
    Mesh N = read(N_filename, orientation);

    for (auto c : N.cells()) {
        if (N.degenerate_or_inverted(c)) {
            std::cout << "TetGen" << std::endl;
            if (!retetrahedrize(N)) {
                return 0;
            }
            break;
        }
    }

    std::vector<Mesh> components = sd(N);
    if (components.size() == 0) {
        return 0;
    }

    /*
    if (out_dir != "") {
        for (int i = 0; i < components.size(); i++) {
            for (int j = 0; j < 2; j++) {
                Mesh& mesh = j == 0 ? components[i].first : components[i].second;
                auto Q = mesh.property<Vertex, Vector3q>("Q");
                auto Q_string = mesh.property<Vertex, std::string>("Q_string");
                for (auto v : mesh.vertices()) {
                    std::stringstream ss;
                    ss << Q[v].transpose();
                    Q_string[v] = ss.str();
                }
                std::stringstream ss;
                ss << out_dir << SEP << (j == 0 ? "M" : "N") << i << ".ovm";
                OpenVolumeMesh::IO::FileManager fm;
                fm.setVerbosityLevel(0);
                fm.writeFile(ss.str(), mesh);
            }
        }
    }
    */

#ifdef GUI
    auto clr = viewer_mesh.property<Cell, Eigen::Vector3d>("clr", {1, 1, 1});
    std::random_device rd;
    std::default_random_engine gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    for (int i = 0; i < components.size(); i++) {
        Mesh& Mi = components[i];
        Eigen::Vector3d color = normalize_color({dis(gen), dis(gen), dis(gen)});
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
    }
    viewer.start(callback, path, "Star Decomposition Maps");
#endif
}
