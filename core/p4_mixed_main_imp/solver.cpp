#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <functional>

#include "solver.h"
#include "../shared/tridiag.h"
#include "../shared/grid.h"

static const double XI_JUMP = M_PI / 4.0;

// Краевые условия 3-го рода: -k*u'(0) + γ₁*u(0) = θ₁,  k*u'(1) + γ₂*u(1) = θ₂
static const double GAMMA1  = 1.0;
static const double THETA1  = 2.0;
static const double GAMMA2  = 3.0;
static const double THETA2  = 4.0;
static const double MU1     = 1.0;   // Теперь это значение u(0) - Дирихле
static const double MU2     = 0.0;   // Теперь это поток w(1) - Нейман

// Функции коэффициентов (Вариант 5)
double k1(double x) { return std::sqrt(2.0) * std::sin(x) + 0.1; }
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

    // 1. Левое ГУ (i = 0): 3-го рода, -k*u'(0) + γ₁*u(0) = θ₁
    //    Аппроксимация балансом на [0, h/2]:
    //    (k_{1/2}/h)*(u₀ - u₁) + (γ₁ + q̄₀*(h/2))*u₀ = f̄₀*(h/2) + θ₁
    double k_half_left = grid::k_half(k1, 0.0, g.half(0));
    double q_bar_0     = grid::q_bar(q1, 0.0, g.half(0));
    double f_bar_0     = grid::f_bar(f1, 0.0, g.half(0));

    a[0] = 0.0;
    b[0] = (k_half_left / h) + GAMMA1 + q_bar_0 * (h / 2.0);
    c[0] = -(k_half_left / h);
    d[0] = f_bar_0 * (h / 2.0) + THETA1;
    // 1. Левое ГУ (i = 0): Дирихле u(0) = MU1
    // Аппроксимация: 1*u_0 + 0*u_1 = MU1
    a[0] = 0.0;
    b[0] = 1.0;
    c[0] = 0.0;
    d[0] = MU1;

    int m = grid::jump_node(g, XI_JUMP);

    // 2. Внутренние узлы (i = 1 .. n-1)
    for (int i = 1; i < n; ++i) {
        double k_left, k_right, q_cell, f_cell;

        // k_left (между x_{i-1} и x_i)
        if (g.x[i] <= XI_JUMP) 
            k_left = grid::k_half(k1, g.x[i-1], g.x[i]);
        else if (g.x[i-1] >= XI_JUMP) 
            k_left = grid::k_half(k2, g.x[i-1], g.x[i]);
        else 
            k_left = grid::k_half_jump(k1, k2, g.x[i-1], g.x[i], XI_JUMP);

        // k_right (между x_i и x_{i+1})
        if (g.x[i+1] <= XI_JUMP) 
            k_right = grid::k_half(k1, g.x[i], g.x[i+1]);
        else if (g.x[i] >= XI_JUMP) 
            k_right = grid::k_half(k2, g.x[i], g.x[i+1]);
        else 
            k_right = grid::k_half_jump(k1, k2, g.x[i], g.x[i+1], XI_JUMP);

        // Интегральные q и f
        if (i < m) {
            q_cell = grid::q_bar(q1, g.left_half(i), g.right_half(i));
            f_cell = grid::f_bar(f1, g.left_half(i), g.right_half(i));
        } else if (i > m) {
            q_cell = grid::q_bar(q2, g.left_half(i), g.right_half(i));
            f_cell = grid::f_bar(f2, g.left_half(i), g.right_half(i));
        } else {
            q_cell = grid::q_bar_jump(q1, q2, g.left_half(i), g.right_half(i), XI_JUMP);
            f_cell = grid::f_bar_jump(f1, f2, g.left_half(i), g.right_half(i), XI_JUMP);
        }

        a[i] = -k_left / h;
        c[i] = -k_right / h;
        b[i] = -(a[i] + c[i]) + q_cell * h;
        d[i] = f_cell * h;
    }

    // 3. Правое ГУ (i = n): 3-го рода, k*u'(1) + γ₂*u(1) = θ₂
    //    Аппроксимация балансом на [x_{n-1/2}, 1]:
    //    (k_{n-1/2}/h)*(uₙ - u_{n-1}) + (γ₂ + q̄ₙ*(h/2))*uₙ = f̄ₙ*(h/2) + θ₂
    double k_half_right = grid::k_half(k2, g.x[n - 1], g.x[n]);
    double q_bar_n      = grid::q_bar(q2, g.half(n - 1), 1.0);
    double f_bar_n      = grid::f_bar(f2, g.half(n - 1), 1.0);

    a[n] = -(k_half_right / h);
    b[n] = (k_half_right / h) + GAMMA2 + q_bar_n * (h / 2.0);
    c[n] = 0.0;
    d[n] = f_bar_n * (h / 2.0) + THETA2;

    return BVP_Result(g, tridiag::solve(a, b, c, d),
                      XI_JUMP, GAMMA1, THETA1, GAMMA2, THETA2,
                      "robin", "robin");
    // 3. Правое ГУ (i = n): Улучшенный Нейман w(1) = MU2
    // Условие: -k(1)*u'(1) = MU2. Интегрируем баланс на [x_{n-1/2}, 1]
    double k_last = grid::k_half(k2, g.left_half(n), 1.0);
    double q_last = grid::q_bar(q2, g.left_half(n), 1.0);
    double f_last = grid::f_bar(f2, g.left_half(n), 1.0);

    // Уравнение: -(k_last/h)*u_{n-1} + (k_last/h + q_last*h/2)*u_n = f_last*h/2 - MU2
    a[n] = -k_last / h;
    b[n] = (k_last / h) + q_last * (h / 2.0);
    c[n] = 0.0;
    d[n] = f_last * (h / 2.0) - MU2;

    return BVP_Result(g, tridiag::solve(a, b, c, d),
                      XI_JUMP, MU1, MU2,
                      "dirichlet", "neumann");
}