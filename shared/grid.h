#pragma once
#include <vector>
#include <functional>
#include <stdexcept>
#include <cmath>

/**
 * grid.h — равномерная сетка и вспомогательные интегралы
 *          для метода баланса (интегрально-интерполяционный метод)
 *
 * Сетка:    x[i]     = i * h,  h = 1/n,  i = 0 .. n
 * Полуузлы: x[i+1/2] = (i + 0.5) * h
 *
 * Используется всеми членами команды — подключение только через #include.
 */

namespace grid {

// ─────────────────────────────────────────────
// 1. Структура равномерной сетки
// ─────────────────────────────────────────────

struct UniformGrid {
    int    n; // число разбиений (узлов n + 1) 
    double h; // шаг сетки h = 1 / n
    std::vector<double> x; // узлы x[i] = i * h, i = 0..n

    explicit UniformGrid(int n_) : n(n_), h(1.0 / n_), x(n_ + 1) {
        if (n_ < 2)
            throw std::runtime_error("grid::UniformGrid: n must be >= 2");
        for (int i = 0; i <= n_; ++i)
            x[i] = i * h;
    } 

    // Полуузел x[i + 0.5]
    double half(int i) const {
        return (i + 0.5) * h;
    }

    // Левый  край ячейки узла i: x[i - 0.5]  (i >= 1)
    double left_half(int i) const {
        return (i - 0.5) * h;
    }

    // Правый край ячейки узла i: x[i + 0.5]  (i <= n-1)
    double right_half(int i) const {
        return (i + 0.5) * h;
    }
};

// ─────────────────────────────────────────────
// 2. Численное интегрирование (правило Симпсона)
// ─────────────────────────────────────────────
 
/**
 * integrate_simpson — интеграл функции f на [a, b] методом Симпсона.
 * steps — число разбиений (чётное), по умолчанию 100.
 * Точность O(h^4) — достаточно для коэффициентов схемы 2-го порядка.
 */
inline double integrate_simpson(
    const std::function<double(double)>& f,
    double a, double b,
    int steps = 100)
{
    if (steps % 2 != 0) ++steps;
    double dt = (b - a) / steps;
    double s = f(a) + f(b);
    for (int i = 1; i < steps; ++i)
        s += f(a + i * dt) * (i % 2 == 0 ? 2.0 : 4.0);
    return s * dt / 3.0;
}
 
// ─────────────────────────────────────────────
// 3. Коэффициенты схемы метода баланса
//    для «обычных» узлов (ячейка не пересекает ξ)
// ─────────────────────────────────────────────
 
/**
 * Интегральное среднее k на полуинтервале [x_i, x_{i+1}]:
 *   k_hat_{i+1/2} = (1/h) * ∫_{x_i}^{x_{i+1}} k(x) dx
 *
 * Используется как коэффициент a_{i+1} = k_hat_{i+1/2} / h²
 * (или c_i = k_hat_{i+1/2} / h²  в зависимости от обозначений).
 */

inline double k_half(
    const std::function<double(double)>& k,
    double xi, double xi1,
    int steps = 100)
{
    double h = xi1 - xi;
    // ∫(1/k)dx по [xi, xi1]
    auto inv_k = [&k](double x){ return 1.0 / k(x); };
    double integral = integrate_simpson(inv_k, xi, xi1, steps);
    return h / integral;  // h / ∫(1/k)dx
}
 
/**
 * Интегральное среднее q на ячейке [x_{i-1/2}, x_{i+1/2}]:
 *   q_bar_i = (1/h) * ∫_{x_{i-1/2}}^{x_{i+1/2}} q(x) dx
 */
inline double q_bar(
    const std::function<double(double)>& q,
    double x_left_half, double x_right_half,
    int steps = 100)
{
    double h = x_right_half - x_left_half;
    return integrate_simpson(q, x_left_half, x_right_half, steps) / h;
}
 
/**
 * Интегральное среднее f на ячейке [x_{i-1/2}, x_{i+1/2}]:
 *   f_bar_i = (1/h) * ∫_{x_{i-1/2}}^{x_{i+1/2}} f(x) dx
 */
inline double f_bar(
    const std::function<double(double)>& f,
    double x_left_half, double x_right_half,
    int steps = 100)
{
    double h = x_right_half - x_left_half;
    return integrate_simpson(f, x_left_half, x_right_half, steps) / h;
}
 
// ─────────────────────────────────────────────
// 4. Коэффициенты схемы в окрестности точки разрыва ξ
//    Ячейка узла m содержит ξ: x_{m-1/2} < ξ < x_{m+1/2}
// ─────────────────────────────────────────────
 
/**
 * k_half_jump — k_hat на полуинтервале [x_i, x_{i+1}],
 * когда ξ ∈ (x_i, x_{i+1}):
 *   k_hat = (1/h) * [ ∫_{x_i}^{ξ} k1(x)dx + ∫_{ξ}^{x_{i+1}} k2(x)dx ]
 *
 * Используется для вычисления a_m и c_{m-1} вблизи разрыва.
 */

inline double k_half_jump(
    const std::function<double(double)>& k1,
    const std::function<double(double)>& k2,
    double xi_node, double xi1_node,
    double xi_jump,
    int steps = 100)
{
    double h = xi1_node - xi_node;
    auto inv_k1 = [&k1](double x){ return 1.0 / k1(x); };
    auto inv_k2 = [&k2](double x){ return 1.0 / k2(x); };
    double left  = integrate_simpson(inv_k1, xi_node,  xi_jump,   steps);
    double right = integrate_simpson(inv_k2, xi_jump,  xi1_node,  steps);
    return h / (left + right);
}
 
/**
 * q_bar_jump — q_bar для особой ячейки узла m,
 * когда ξ ∈ (x_{m-1/2}, x_{m+1/2}):
 *   q_bar_m = (1/h) * [ ∫_{x_{m-1/2}}^{ξ} q1(x)dx + ∫_{ξ}^{x_{m+1/2}} q2(x)dx ]
 */
inline double q_bar_jump(
    const std::function<double(double)>& q1,
    const std::function<double(double)>& q2,
    double x_left_half, double x_right_half,
    double xi_jump,
    int steps = 100)
{
    double h = x_right_half - x_left_half;
    double left_part  = integrate_simpson(q1, x_left_half, xi_jump,      steps);
    double right_part = integrate_simpson(q2, xi_jump,     x_right_half, steps);
    return (left_part + right_part) / h;
}
 
/**
 * f_bar_jump — f_bar для особой ячейки узла m,
 * когда ξ ∈ (x_{m-1/2}, x_{m+1/2}):
 *   f_bar_m = (1/h) * [ ∫_{x_{m-1/2}}^{ξ} f1(x)dx + ∫_{ξ}^{x_{m+1/2}} f2(x)dx ]
 */
inline double f_bar_jump(
    const std::function<double(double)>& f1,
    const std::function<double(double)>& f2,
    double x_left_half, double x_right_half,
    double xi_jump,
    int steps = 100)
{
    double h = x_right_half - x_left_half;
    double left_part  = integrate_simpson(f1, x_left_half, xi_jump,      steps);
    double right_part = integrate_simpson(f2, xi_jump,     x_right_half, steps);
    return (left_part + right_part) / h;
}
 
// ─────────────────────────────────────────────
// 5. Вспомогательное: найти индекс особого узла m
//    (узел, чья ячейка содержит точку разрыва ξ)
// ─────────────────────────────────────────────
 
/**
 * jump_node — возвращает индекс m такой, что
 *   x[m-1/2] < xi_jump < x[m+1/2]
 * то есть m = round(xi_jump / h), ограниченный [1, n-1].
 *
 * Если ξ точно совпадает с узлом сетки — возвращает этот узел
 * (интегралы корректно считаются как предельный случай).
 */
inline int jump_node(const UniformGrid& g, double xi_jump) {
    int m = static_cast<int>(std::round(xi_jump / g.h));
    if (m < 1)   m = 1;
    if (m > g.n - 1) m = g.n - 1;
    return m;
}
 
// ─────────────────────────────────────────────
// 6. Нормы сеточных функций
// ─────────────────────────────────────────────
 
/**
 * max_norm — максимум |u[i] - v[i]| по всем узлам i=0..n.
 * Используется для вычисления ε1 и ε2.
 *
 * Возвращает {max_value, index_of_max}.
 */
inline std::pair<double, int> max_norm(
    const std::vector<double>& u,
    const std::vector<double>& v)
{
    if (u.size() != v.size())
        throw std::runtime_error("grid::max_norm: векторы разной длины");
    double mx = 0.0;
    int idx = 0;
    for (int i = 0; i < static_cast<int>(u.size()); ++i) {
        double diff = std::abs(u[i] - v[i]);
        if (diff > mx) { mx = diff; idx = i; }
    }
    return {mx, idx};
}
 
/**
 * max_norm_coarse — максимум |v1[i] - v2[2*i]| по общим узлам
 * сетки n и сетки 2n. Используется для вычисления ε2.
 *
 * v1 — решение на сетке n  (размер n+1)
 * v2 — решение на сетке 2n (размер 2n+1)
 *
 * Возвращает {max_value, index_of_max} (индекс в сетке n).
 */
inline std::pair<double, int> max_norm_coarse(
    const std::vector<double>& v1,
    const std::vector<double>& v2)
{
    int n = static_cast<int>(v1.size()) - 1;
    if (static_cast<int>(v2.size()) != 2 * n + 1)
        throw std::runtime_error("grid::max_norm_coarse: v2 должен быть размера 2n+1");
    double mx = 0.0;
    int idx = 0;
    for (int i = 0; i <= n; ++i) {
        double diff = std::abs(v1[i] - v2[2 * i]);
        if (diff > mx) { mx = diff; idx = i; }
    }
    return {mx, idx};
}

} // namespace grid
