#pragma once
#include <vector>
#include <string>
#include "../shared/grid.h"

struct BVP_Result {
    grid::UniformGrid grid;
    std::vector<double> u;

    // Параметры задачи для фронтенда
    double xi_jump;
    double gamma1;
    double theta1;
    double gamma2;
    double theta2;
    std::string bc_left;
    std::string bc_right;

    BVP_Result() : grid(2), u(), xi_jump(0.0),
                   gamma1(0.0), theta1(0.0), gamma2(0.0), theta2(0.0),
                   bc_left("unknown"), bc_right("unknown") {}

    BVP_Result(const grid::UniformGrid& g, const std::vector<double>& u_vec,
               double xi, double g1, double t1, double g2, double t2,
               const std::string& bcl, const std::string& bcr)
        : grid(g), u(u_vec), xi_jump(xi),
          gamma1(g1), theta1(t1), gamma2(g2), theta2(t2),
          bc_left(bcl), bc_right(bcr) {}
};

BVP_Result solve_mixed_main_imp(int n);