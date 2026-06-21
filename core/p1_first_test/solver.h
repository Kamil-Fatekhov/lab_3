#pragma once
#include <vector>
#include <cmath>
#include "math_utils.h"

struct TaskConfig {
    double xi = M_PI / 4.0;
    double mu1 = 1.0;
    double mu2 = 0.0;
    double k1 = 1.0;
    double k2 = 0.5;
    double q1 = 1.0;
    double q2 = (M_PI * M_PI) / 16.0;
    double f1 = 1.0;
    double f2 = M_SQRT2 / 2.0;
    double omega;
    std::vector<double> C; 

    TaskConfig() {
        omega = std::sqrt(q2 / k2);
        std::vector<std::vector<double>> A = {
            {1.0, 1.0, 0.0, 0.0},
            {0.0, 0.0, std::exp(omega), std::exp(-omega)},
            {std::exp(xi), std::exp(-xi), -std::exp(omega * xi), -std::exp(-omega * xi)},
            {k1 * std::exp(xi), -k1 * std::exp(-xi), -k2 * omega * std::exp(omega * xi), k2 * omega * std::exp(-omega * xi)}
        };
        std::vector<double> B = {
            0.0,
            -f2 / q2,
            f2 / q2 - 1.0,
            0.0
        };
        C = solve_4x4(A, B);
    }

    double exact_solution(double x) const {
        if (x < xi) {
            return C[0] * std::exp(x) + C[1] * std::exp(-x) + (f1 / q1);
        } else {
            return C[2] * std::exp(omega * x) + C[3] * std::exp(-omega * x) + (f2 / q2);
        }
    }
};

class BalanceSolver {
private:
    TaskConfig cfg;
    int n;
    double h;

    double calc_a(int i) {
        double x_prev = (i - 1) * h;
        double x_curr = i * h;
        if (cfg.xi > x_prev && cfg.xi <= x_curr) {
            return 1.0 / (((cfg.xi - x_prev) / (h * cfg.k1)) + ((x_curr - cfg.xi) / (h * cfg.k2)));
        } else if (x_curr <= cfg.xi) {
            return cfg.k1;
        } else {
            return cfg.k2;
        }
    }

    double calc_d(int i) {
        double x_left = (i - 0.5) * h;
        double x_right = (i + 0.5) * h;
        if (cfg.xi > x_left && cfg.xi <= x_right) {
            return (cfg.q1 * (cfg.xi - x_left) + cfg.q2 * (x_right - cfg.xi)) / h;
        } else if (x_right <= cfg.xi) {
            return cfg.q1;
        } else {
            return cfg.q2;
        }
    }

    double calc_phi(int i) {
        double x_left = (i - 0.5) * h;
        double x_right = (i + 0.5) * h;
        if (cfg.xi > x_left && cfg.xi <= x_right) {
            return (cfg.f1 * (cfg.xi - x_left) + cfg.f2 * (x_right - cfg.xi)) / h;
        } else if (x_right <= cfg.xi) {
            return cfg.f1;
        } else {
            return cfg.f2;
        }
    }

public:
    BalanceSolver(int grids) : n(grids), h(1.0 / grids) {}

    std::vector<double> solve() {
        std::vector<double> a(n + 1), c(n + 1), b(n + 1), f(n + 1);
        for (int i = 1; i < n; ++i) {
            double ai = calc_a(i);
            double ai_plus_1 = calc_a(i + 1);
            a[i] = ai / (h * h);
            b[i] = ai_plus_1 / (h * h);
            c[i] = a[i] + b[i] + calc_d(i);
            f[i] = calc_phi(i);
        }

        std::vector<double> alpha(n + 1), beta(n + 1);
        alpha[1] = 0.0;
        beta[1] = cfg.mu1;

        for (int i = 1; i < n; ++i) {
            double denom = c[i] - a[i] * alpha[i];
            alpha[i + 1] = b[i] / denom;
            beta[i + 1] = (f[i] + a[i] * beta[i]) / denom;
        }

        std::vector<double> v(n + 1);
        v[n] = cfg.mu2;
        for (int i = n - 1; i >= 0; --i) {
            v[i] = alpha[i + 1] * v[i + 1] + beta[i + 1];
        }
        return v;
    }

    double max_error(const std::vector<double>& v) {
        double max_err = 0.0;
        for (int i = 0; i <= n; ++i) {
            double x = i * h;
            double u_exact = cfg.exact_solution(x);
            max_err = std::max(max_err, std::abs(u_exact - v[i]));
        }
        return max_err;
    }
};

// В самый конец файла core/p1_first_test/solver.h

inline std::vector<double> solve_p1_task(int n) {
    BalanceSolver solver(n);
    return solver.solve();
}

inline double get_analytical_p1(double x) {
    TaskConfig cfg;
    return cfg.exact_solution(x);
}