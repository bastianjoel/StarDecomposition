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

std::pair<StarCenterResult, Eigen::Vector3d> star_center_close_to(const Eigen::Vector3d& point, const std::vector<Eigen::Vector3d>& positions, const std::vector<Eigen::Vector3d>& normals) {
    CGAL::Quadratic_program<double> lp(CGAL::SMALLER, false);
    lp.set_c(3, 1);
    lp.set_c(4, 1);
    lp.set_c(5, 1);
    for (int i = 0; i < positions.size(); i++) {
        Eigen::Vector3d n = normals[i];
        // n_0 * x + n_1 * y + n_2 * z <= -d
        for (int j = 0; j < 3; j++) {
            lp.set_a(j, i, n(j));
        }
        lp.set_b(i, n.dot(positions[i]));
    }
    for (int i = 0; i < 3; i++) {
        // (x, y, z) - point_i <= u_i
        lp.set_a(i, positions.size() + i, 1);
        lp.set_a(i + 3, positions.size() + i, 1);
        lp.set_b(positions.size() + i, point(i));
        // -((x, y, z) - point) <= u_i
        lp.set_a(i, positions.size() + 3 + i, -1);
        lp.set_a(i + 3, positions.size() + 3 + i, -1);
        lp.set_b(positions.size() + 3 + i, -point(i));
    }
    CGAL::Quadratic_program_solution<double> sol = CGAL::solve_linear_program(lp, 21.0);
    if (sol.is_infeasible()) {
        return std::make_pair(INVALID, Eigen::Vector3d::Zero());
    }

    Eigen::Vector3d x;
    auto it = sol.variable_values_begin();
    for (int i = 0; i < 3; i++) {
        x(i) = CGAL::quotient_truncation(*it);
        it++;
    }

    std::cout << "Vars: " << CGAL::quotient_truncation(*it) << " " << CGAL::quotient_truncation(*(it + 1)) << " " << CGAL::quotient_truncation(*(it + 2));
    std::cout << " " << CGAL::quotient_truncation(*(it + 3)) << " " << CGAL::quotient_truncation(*(it + 4)) << " " << CGAL::quotient_truncation(*(it + 5)) << std::endl;
    std::cout << "Objective: " << sol.objective_value() << std::endl;
    // std::cout << "Solution: " << sol << std::endl;
    std::cout << "Solves: " << sol.is_valid() << std::endl;
    if (sol.objective_value() == 0) {
        return std::make_pair(VALID_EQUAL, point);
    }
    return std::make_pair(VALID, x);
}
