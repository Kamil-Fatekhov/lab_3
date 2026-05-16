#pragma once
#include <vector>
#include <stdexcept>
#include <cmath>
#include <string>

namespace tridiag {

// Система: a[i]*x[i-1] + b[i]*x[i] + c[i]*x[i+1] = d[i]
// Граничные строки задаются снаружи через b[0],d[0] и b[n-1],d[n-1]
inline std::vector<double> solve(
    const std::vector<double>& a,  // нижняя диагональ (a[0] не используется)
    const std::vector<double>& b,  // главная диагональ
    const std::vector<double>& c,  // верхняя диагональ (c[n-1] не используется)
    const std::vector<double>& d)  // правая часть
{
    int n = (int)b.size();
    if (n == 0) return {};

    std::vector<double> alpha(n), beta(n), x(n);

    alpha[0] = -c[0] / b[0];
    beta[0]  =  d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double denom = b[i] + a[i] * alpha[i - 1];
        if (std::abs(denom) < 1e-15)
            throw std::runtime_error("tridiag::solve: нулевой пивот на шаге "
                                     + std::to_string(i));
        alpha[i] = (i < n - 1) ? -c[i] / denom : 0.0;
        beta[i]  = (d[i] - a[i] * beta[i - 1]) / denom;
    }

    x[n - 1] = beta[n - 1];
    for (int i = n - 2; i >= 0; --i)
        x[i] = alpha[i] * x[i + 1] + beta[i];

    return x;
}

// Обёртка для метода баланса: -A[i]*v[i-1] + B[i]*v[i] - C[i]*v[i+1] = D[i]
// mu1 = v[0], mu2 = v[n-1] — граничные условия первого рода
inline std::vector<double> solve_balance(
    const std::vector<double>& A,
    const std::vector<double>& B,
    const std::vector<double>& C,
    const std::vector<double>& D,
    double mu1, double mu2)
{
    int n = (int)B.size();
    std::vector<double> a(n,0), b(n), c(n,0), d(n);

    b[0] = 1.0; d[0] = mu1;

    for (int i = 1; i < n - 1; ++i) {
        a[i] = -A[i];
        b[i] =  B[i];
        c[i] = -C[i];
        d[i] =  D[i];
    }

    b[n-1] = 1.0; d[n-1] = mu2;

    return solve(a, b, c, d);
}

} // namespace tridiag