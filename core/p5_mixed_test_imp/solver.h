#pragma once
#include <vector>

// Функция вычисления точного аналитического решения в точке x
double get_analytical_test(double x);

// Функция численного решения смешанной краевой задачи на сетке из n интервалов
// Реализует метод баланса с улучшенной аппроксимацией ГУ третьего рода O(h^2)
std::vector<double> solve_test_mixed(int n);