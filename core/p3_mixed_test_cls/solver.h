#pragma once
#include <vector>

// Численное решение тестовой задачи (классическая аппроксимация 3-го рода)
std::vector<double> solve_test_mixed(int n);

// Точное аналитическое решение в точке x
double get_analytical_test(double x);