#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include "solver.h"
#include "grid.h"

int main() {
    // Настройка локали для корректного вывода символов
    std::setlocale(LC_ALL, "Russian");

    std::cout << "=== BVP2 Mixed Test Improved Approximation O(h^2) (Variant 5) ===\n";
    std::cout << "Режим подробного вывода (печать каждого узла сетки)\n\n";
    
    int n = 0;
    std::cout << "Введите число разбиений сетки (n): ";
    std::cin >> n;

    // Валидация ввода
    if (n <= 1) {
        std::cerr << "Ошибка: число разбиений n должно быть больше 1!\n";
        return 1;
    }
    
    // 1. Находим численное решение v на заданной сетке
    std::vector<double> v = solve_test_mixed(n);
    grid::UniformGrid g(n);
    
    // 2. Считаем точное аналитическое решение u
    std::vector<double> u(n + 1);
    for (int i = 0; i <= n; ++i) {
        u[i] = get_analytical_test(g.x[i]);
    }
    
    // 3. Находим максимальное абсолютное отклонение и его индекс
    auto norm = grid::max_norm(u, v);
    double eps1 = norm.first;
    int max_err_idx = norm.second;

    // Вывод краткой справки
    std::cout << "\nСправка по решению тестовой задачи для n = " << n << ":\n";
    std::cout << "-------------------------------------------------------\n";
    std::cout << "Шаг равномерной сетки h                = " << std::fixed << std::setprecision(8) << g.h << "\n";
    std::cout << "Максимальная погрешность eps1          = " << std::scientific << eps1 << "\n";
    std::cout << "Максимальная ошибка наблюдается в узле = " << max_err_idx << " при x = " 
              << std::fixed << std::setprecision(6) << g.x[max_err_idx] << "\n";
    std::cout << "-------------------------------------------------------\n\n";

    // Полный вывод таблицы результатов (шаг строго равен 1)
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::left << std::setw(10) << "№ узла" 
              << std::setw(15) << "Координата x_i" 
              << std::setw(18) << "Точное u(x_i)" 
              << std::setw(18) << "Численное v(x_i)" 
              << std::setw(20) << "Разность u_i - v_i" << "\n";
    std::cout << std::string(80, '-') << "\n";
    
    // Цикл идет по абсолютно каждому элементу без пропусков
    for (int i = 0; i <= n; ++i) {
        std::cout << std::left << std::setw(10) << i 
                  << std::setw(15) << std::fixed << std::setprecision(6) << g.x[i]
                  << std::setw(18) << std::setprecision(10) << u[i]
                  << std::setw(18) << std::setprecision(10) << v[i]
                  << std::setw(20) << std::scientific << (u[i] - v[i]) << "\n";
    }
    std::cout << std::string(80, '-') << "\n";

    return 0;
}