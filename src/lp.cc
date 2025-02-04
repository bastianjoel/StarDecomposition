#include "lp.h"
#include <utility>

Eigen::Vector3d kernel_chebyshev_center(const std::vector<Eigen::Vector3d>& positions, const std::vector<Eigen::Vector3d>& normals) {
    CGAL::Quadratic_program<double> lp(CGAL::SMALLER, false);
    lp.set_c(3, -1);
    for (int i = 0; i < positions.size(); i++) {
        Eigen::Vector3d n = normals[i];
        // n_0 * x + n_1 * y + n_2 * z + r <= -d
        for (int j = 0; j < 3; j++) {
            lp.set_a(j, i, n(j));
        }
        lp.set_a(3, i, 1);
        lp.set_b(i, n.dot(positions[i]));
    }
    CGAL::Quadratic_program_solution<double> sol = CGAL::solve_linear_program(lp, 21.0);
    Eigen::Vector3d x;
    auto it = sol.variable_values_begin();
    for (int i = 0; i < 3; i++) {
        x(i) = CGAL::quotient_truncation(*it);
        it++;
    }
    if (sol.is_optimal() && !std::isnan(x(0)) && !std::isnan(x(1)) && !std::isnan(x(2))) {
        return x;
    }
    x = Eigen::Vector3d::Zero();
    for (int i = 0; i < positions.size(); i++) {
        x += positions[i];
    }
    x /= (double) positions.size();
    return x;
}

std::pair<StarCenterResult, Eigen::Vector3d> star_center_close_to(const Eigen::Vector3d& normal, const std::vector<Eigen::Vector3d>& positions, const std::vector<Eigen::Vector3d>& normals) {
    CGAL::Quadratic_program<double> lp(CGAL::SMALLER, false);
    lp.set_c(0, -normal[0]);
    lp.set_c(1, -normal[1]);
    lp.set_c(2, -normal[2]);
    for (int i = 0; i < positions.size(); i++) {
        Eigen::Vector3d n = normals[i];
        // n_0 * x + n_1 * y + n_2 * z <= -d
        for (int j = 0; j < 3; j++) {
            lp.set_a(j, i, n(j));
        }
        lp.set_b(i, n.dot(positions[i]));
    }

    CGAL::Quadratic_program_solution<double> sol = CGAL::solve_linear_program(lp, 21.0);
    Eigen::Vector3d x;
    auto it = sol.variable_values_begin();
    for (int i = 0; i < 3; i++) {
        x(i) = CGAL::quotient_truncation(*it);
        it++;
    }

    it = sol.variable_values_begin();
    std::cout << "Vars: " << CGAL::quotient_truncation(*it) << " " << CGAL::quotient_truncation(*(it + 1)) << " " << CGAL::quotient_truncation(*(it + 2)) << std::endl;
    std::cout << "Objective: " << sol.objective_value() << std::endl;
    std::cout << "Point: " << normal[0] << " " << normal[1] << " " << normal[2] << std::endl;
    // std::cout << "Solution: " << sol << std::endl;
    std::cout << "Solves: " << sol.is_valid() << std::endl;
    if (sol.is_infeasible()) {
        std::cout << "Infeasible" << std::endl;
        return std::make_pair(INVALID, Eigen::Vector3d::Zero());
    }

    return std::make_pair(VALID, x);
}
