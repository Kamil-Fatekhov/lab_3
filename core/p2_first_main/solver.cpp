#include "solver.h"
#include "../../shared/tridiag.h"
#include "../../shared/grid.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <algorithm>

// ── Параметры задачи ─────────────────────────────────────────────────────────

static const double XI = M_PI / 4.0; // точка разрыва
static const double MU1 = 1.0;       // u(0) = mu1
static const double MU2 = 0.0;       // u(1) = mu2

// ── Коэффициенты ─────────────────────────────────────────────────────────────

double k1(double x) { return std::sqrt(2.0) * std::sin(x); }
double k2(double x) { return std::cos(x) * std::cos(x); }
double q1(double x) { (void)x; return 1.0; }
double q2(double x) { return x * x; }
double f1(double x) { return std::sin(2.0 * x); }
double f2(double x) { return std::cos(x); }

double k(double x) { return x < XI ? k1(x) : k2(x); }
double q(double x) { return x < XI ? q1(x) : q2(x); }
double f(double x) { return x < XI ? f1(x) : f2(x); }

// ── Построение и решение схемы ───────────────────────────────────────────────

/*
std::vector<double> solve_bvp(int n) {
    grid::UniformGrid g(n);
    const double h = g.h; // шаг сетки

    // Находим индекс особого узла m: ячейка [x_{m-1/2}, x_{m+1/2}] содержит XI
    int m = grid::jump_node(g, XI);

    // Коэффициенты трёхдиагональной системы схемы метода баланса:
    //   -A[i]*v[i-1] + B[i]*v[i] - C[i]*v[i+1] = D[i]
    std::vector<double> A(n + 1, 0.0);
    std::vector<double> B(n + 1, 0.0);
    std::vector<double> C(n + 1, 0.0);
    std::vector<double> D(n + 1, 0.0);

    for (int i = 1; i < n; ++i) {
        double xl = g.left_half(i); // x_{i-1/2}
        double xr = g.right_half(i); // x_{i+1/2}

        double k_left, k_right, q_bar_i, f_bar_i;

        if (i == m) {
            // Особый узел: ячейки и полуинтервалы пересекают точку разрыва XI
 
            // k_hat_{i-1/2} = (1/h) * [∫_{x_{i-1}}^{XI} k1 dx + ∫_{XI}^{x_i} k2 dx]
            // Но полуинтервал для A[i] — это [x_{i-1}, x_i], где x_{i-1} = (i-1)*h
            // Аналогично для C[i] — полуинтервал [x_i, x_{i+1}]
            // Проверяем в каком полуинтервале лежит XI

            double x_im1 = g.x[i - 1];  // x_{i-1}
            double x_i   = g.x[i];      // x_i
            double x_ip1 = g.x[i + 1];  // x_{i+1}

            // Левый полуинтервал [x[i-1], x[i]] — XI может быть внутри или нет
            if (XI > x_im1 && XI < x_i) {
                k_left = grid::k_half_jump(k1, k2, x_im1, x_i, XI);
            } else {
                double mid = (x_im1 + x_i) * 0.5;
                k_left = grid::k_half(mid < XI ? k1 : k2, x_im1, x_i);
            }

            // Правый полуинтервал [x[i], x[i+1]]
            if (XI > x_i && XI < x_ip1) {
                k_right = grid::k_half_jump(k1, k2, x_i, x_ip1, XI);
            } else {
                double mid = (x_i + x_ip1) * 0.5;
                k_right = grid::k_half(mid < XI ? k1 : k2, x_i, x_ip1);
            }

            // Ячейка [xl, xr] — XI точно внутри (по определению jump_node)
            q_bar_i = grid::q_bar_jump(q1, q2, xl, xr, XI);
            f_bar_i = grid::f_bar_jump(f1, f2, xl, xr, XI);

        } else {
            // Левый полуинтервал [x[i-1], x[i]]
            // Правый полуинтервал [x[i], x[i+1]]
            // Выбираем ветвь по середине каждого полуинтервала

            double mid_left  = (g.x[i-1] + g.x[i])     * 0.5;
            double mid_right = (g.x[i]   + g.x[i+1])   * 0.5;

            auto k_left_func  = (mid_left  < XI) ? k1 : k2;
            auto k_right_func = (mid_right < XI) ? k1 : k2;
            auto q_func       = (g.x[i]   < XI) ? q1 : q2;
            auto f_func       = (g.x[i]   < XI) ? f1 : f2;

            k_left  = grid::k_half(k_left_func,  g.x[i-1], g.x[i]);
            k_right = grid::k_half(k_right_func, g.x[i],   g.x[i+1]);
            q_bar_i = grid::q_bar(q_func, xl, xr);
            f_bar_i = grid::f_bar(f_func, xl, xr);   
        }

        A[i] = k_left  / h;
        C[i] = k_right / h;
        B[i] = A[i] + C[i] + q_bar_i * h;
        D[i] = f_bar_i * h;

    }

    return tridiag::solve_balance(A, B, C, D, MU1, MU2);
}
*/
std::vector<double> solve_bvp(int n) {
    grid::UniformGrid g(n);
    const double h = g.h;
    int m = grid::jump_node(g, XI);

    std::vector<double> A(n+1, 0.0);
    std::vector<double> B(n+1, 0.0);
    std::vector<double> C(n+1, 0.0);
    std::vector<double> D(n+1, 0.0);

    for (int i = 1; i < n; ++i) {
        double k_left, k_right, q_cell, f_cell;

        // k_left — на полуинтервале [x[i-1], x[i]]
        if (g.x[i] <= XI)
            k_left = grid::k_half(k1, g.x[i-1], g.x[i]);
        else if (g.x[i-1] >= XI)
            k_left = grid::k_half(k2, g.x[i-1], g.x[i]);
        else
            k_left = grid::k_half_jump(k1, k2, g.x[i-1], g.x[i], XI);

        // k_right — на полуинтервале [x[i], x[i+1]]
        if (g.x[i+1] <= XI)
            k_right = grid::k_half(k1, g.x[i], g.x[i+1]);
        else if (g.x[i] >= XI)
            k_right = grid::k_half(k2, g.x[i], g.x[i+1]);
        else
            k_right = grid::k_half_jump(k1, k2, g.x[i], g.x[i+1], XI);

        // q и f — по ячейке [x[i-1/2], x[i+1/2]]
        if (i < m) {
            q_cell = grid::q_bar(q1, g.left_half(i), g.right_half(i));
            f_cell = grid::f_bar(f1, g.left_half(i), g.right_half(i));
        } else if (i > m) {
            q_cell = grid::q_bar(q2, g.left_half(i), g.right_half(i));
            f_cell = grid::f_bar(f2, g.left_half(i), g.right_half(i));
        } else {
            q_cell = grid::q_bar_jump(q1, q2, g.left_half(i), g.right_half(i), XI);
            f_cell = grid::f_bar_jump(f1, f2, g.left_half(i), g.right_half(i), XI);
        }

        A[i] = k_left  / h;
        C[i] = k_right / h;
        B[i] = A[i] + C[i] + q_cell * h;
        D[i] = f_cell * h;
    }

    return tridiag::solve_balance(A, B, C, D, MU1, MU2);
}

// ── Вычисление ε₂ ────────────────────────────────────────────────────────────

std::pair<double, int> compute_eps2(const std::vector<double>& v, const std::vector<double>& v2) {
    return grid::max_norm_coarse(v, v2);
}


// ── Решение с заданной точностью ─────────────────────────────────────────────

std::pair<std::vector<double>, double> solve_with_accuracy(double eps, int n_start, int& n_out) {
    int n = n_start;
    std::vector<double> v;
    double eps2 = 0.0;
    double prev_eps2 = 1.0e10;

    while (n <= 1000000) {
        v = solve_bvp(n);
        auto v2 = solve_bvp(2 * n);
        auto res = compute_eps2(v, v2);
        eps2 = res.first;

        std::cout << "Debug: n=" << n << " | eps2=" << eps2 << std::endl;

        if (eps2 <= eps) {
            n_out = n;
            return {v, eps2};
        }
        
        if (n >= 100000 || eps2 > prev_eps2) {
            std::cout << "Debug: Stopping due to limit or stagnation." << std::endl;
            break; 
        } 
        
        prev_eps2 = eps2;
        n *= 2;
    }
    
    // Если дошли сюда, значит точность не достигнута, 
    // но мы возвращаем результат на последнем посчитанном n
    n_out = n;
    return {v, eps2};
}