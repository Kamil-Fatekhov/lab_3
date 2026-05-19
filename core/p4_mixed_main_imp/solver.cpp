#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <functional>

#include "solver.h"
#include "../shared/tridiag.h"
#include "../shared/grid.h"

// Параметры 5 варианта — единственное место определения
static const double XI_JUMP = M_PI / 4.0;
static const double MU1     = 1.0;   // Левое ГУ (поток, Нейман)
static const double MU2     = 0.0;   // Правое ГУ (значение, Дирихле)

double k1(double x) { return std::sqrt(2.0) * std::sin(x); }
double k2(double x) { return std::pow(std::cos(x), 2); }

double q1(double x) { return 1.0; }
double q2(double x) { return x * x; }

double f1(double x) { return std::sin(2.0 * x); }
double f2(double x) { return std::cos(x); }

BVP_Result solve_mixed_main_imp(int n) {
    grid::UniformGrid g(n);
    double h = g.h;

    int num_nodes = n + 1;

    std::vector<double> a(num_nodes, 0.0);
    std::vector<double> b(num_nodes, 0.0);
    std::vector<double> c(num_nodes, 0.0);
    std::vector<double> d(num_nodes, 0.0);

    // 1. Левое ГУ (i = 0): Нейман, улучшенная аппроксимация O(h^2)
    double k_half   = grid::k_half(k1, 0.0, g.half(0));
    double q_bar_0  = grid::q_bar(q1, 0.0, g.half(0));
    double f_bar_0  = grid::f_bar(f1, 0.0, g.half(0));

    a[0] = 0.0;
    b[0] = (k_half / h) + q_bar_0 * (h / 2.0);
    c[0] = -(k_half / h);
    d[0] = f_bar_0 * (h / 2.0) + MU1;

    int m = grid::jump_node(g, XI_JUMP);

    // 2. Внутренние узлы (i = 1 .. n-1)
    for (int i = 1; i < n; ++i) {
        double k_left, k_right, q_cell, f_cell;

        if (g.x[i] <= XI_JUMP) {
            k_left = grid::k_half(k1, g.x[i - 1], g.x[i]);
        } else if (g.x[i - 1] >= XI_JUMP) {
            k_left = grid::k_half(k2, g.x[i - 1], g.x[i]);
        } else {
            k_left = grid::k_half_jump(k1, k2, g.x[i - 1], g.x[i], XI_JUMP);
        }

        if (g.x[i + 1] <= XI_JUMP) {
            k_right = grid::k_half(k1, g.x[i], g.x[i + 1]);
        } else if (g.x[i] >= XI_JUMP) {
            k_right = grid::k_half(k2, g.x[i], g.x[i + 1]);
        } else {
            k_right = grid::k_half_jump(k1, k2, g.x[i], g.x[i + 1], XI_JUMP);
        }

        if (i == m) {
            q_cell = grid::q_bar_jump(q1, q2, g.left_half(i), g.right_half(i), XI_JUMP);
            f_cell = grid::f_bar_jump(f1, f2, g.left_half(i), g.right_half(i), XI_JUMP);
        } else if (i < m) {
            q_cell = grid::q_bar(q1, g.left_half(i), g.right_half(i));
            f_cell = grid::f_bar(f1, g.left_half(i), g.right_half(i));
        } else {
            q_cell = grid::q_bar(q2, g.left_half(i), g.right_half(i));
            f_cell = grid::f_bar(f2, g.left_half(i), g.right_half(i));
        }

        a[i] = -k_left / h;
        c[i] = -k_right / h;
        b[i] = -(a[i] + c[i]) + q_cell * h;
        d[i] = f_cell * h;
    }

    // 3. Правое ГУ (i = n): Дирихле
    a[n] = 0.0;
    b[n] = 1.0;
    c[n] = 0.0;
    d[n] = MU2;

    return BVP_Result(g, tridiag::solve(a, b, c, d),
                      XI_JUMP, MU1, MU2,
                      "neumann", "dirichlet");
}