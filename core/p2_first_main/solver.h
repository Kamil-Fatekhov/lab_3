#pragma once 
#include <vector>
#include <utility>

/**
 * solver.h — первая краевая задача, основная задача, вариант 5
 *
 * -(k(x)*u')' + q(x)*u = f(x),  x in (0,1)
 *  u(0) = mu1 = 1,  u(1) = mu2 = 0
 *
 * k1(x) = 2*sin(x),   k2(x) = cos²(x)     при x > xi
 * q1(x) = 1,          q2(x) = x²           при x > xi
 * f1(x) = sin(2x),    f2(x) = cos(x)       при x > xi
 * xi = pi/4
 */
 
// ── Коэффициенты задачи ──────────────────────────────────────────────────────

double k(double x); // теплопроводность 
double q(double x); // коэффициент при u
double f(double x); // правая часть

// Левые и правые ветви (нужны для интегрирования по частям у узла разрыва)
double k1(double x);
double k2(double x);
double q1(double x);
double q2(double x);
double f1(double x);
double f2(double x);
 
// ── Решатель ─────────────────────────────────────────────────────────────────
 
/**
 * solve_bvp — решает задачу на равномерной сетке размерности n.
 *
 * Возвращает вектор v длины n+1: v[i] ≈ u(x_i), x_i = i/n.
 */
std::vector<double> solve_bvp(int n);

/**
 * solve_with_accuracy — решает задачу с заданной точностью eps.
 *
 * Начинает с n_start, удваивает сетку пока
 *   eps2 = max|v(x_i) - v2(x_{2i})| > eps
 *
 * Возвращает:
 *   .first  — решение на финальной сетке n
 *   .second — достигнутая точность eps2
 *
 * Параметр n_out заполняется финальным n.
 */
std::pair<std::vector<double>, double> solve_with_accuracy(double eps, int n_start, int& n_out);

/**
 * compute_eps2 — вычисляет достигнутую точность между сетками n и 2n.
 * Возвращает {eps2, индекс максимума в сетке n}.
 */
std::pair<double, int> compute_eps2(const std::vector<double>& v, const std::vector<double>& v2);