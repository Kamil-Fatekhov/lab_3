#include <iostream>
#include <iomanip>
#include "solver.h"
#include "grid.h"

int main() {
    std::cout << "=== BVP2 Mixed Test Classical Approximation (Variant 5) ===\n\n";
    
    int n = 10;
    double eps1 = 1.0;
    double target_eps = 0.5e-6; // Требуемая точность из методички
    
    std::vector<double> v;
    std::vector<double> u;
    grid::UniformGrid* final_grid = nullptr;
    int max_err_idx = 0;
    
    // Автоматическое сгущение сеток
    while (true) {
        v = solve_test_mixed(n);
        grid::UniformGrid g(n);
        
        u.resize(n + 1);
        for (int i = 0; i <= n; ++i) {
            u[i] = get_analytical_test(g.x[i]);
        }
        
        auto norm = grid::max_norm(u, v);
        eps1 = norm.first;
        max_err_idx = norm.second;
        
        if (eps1 <= target_eps || n >= 1000000) {
            final_grid = new grid::UniformGrid(n);
            break;
        }
        n *= 2; 
    }

    // Вывод официальной справки по ТЗ
    std::cout << "Справка:\n";
    std::cout << "Для решения задачи использована равномерная сетка с числом разбиений n = " << n << "\n";
    std::cout << "Задача решена с погрешностью eps1 = " << std::scientific << eps1 << "\n";
    std::cout << "Максимальное отклонение аналитического и численного решений наблюдается в точке x = " 
              << std::fixed << std::setprecision(6) << final_grid->x[max_err_idx] << "\n\n";

    // Вывод таблицы (разреженный шаг, чтобы не перегружать консоль)
    int step = std::max(1, n / 20); 
    
    std::cout << std::string(75, '-') << "\n";
    std::cout << std::left << std::setw(10) << "№ узла" 
              << std::setw(15) << "x_i" 
              << std::setw(15) << "u(x_i)" 
              << std::setw(15) << "v(x_i)" 
              << std::setw(20) << "u(x_i) - v(x_i)" << "\n";
    std::cout << std::string(75, '-') << "\n";
    
    for (int i = 0; i <= n; i += step) {
        std::cout << std::left << std::setw(10) << i 
                  << std::setw(15) << std::fixed << std::setprecision(6) << final_grid->x[i]
                  << std::setw(15) << u[i]
                  << std::setw(15) << v[i]
                  << std::setw(20) << std::scientific << (u[i] - v[i]) << "\n";
    }
    
    if (n % step != 0) {
        std::cout << std::left << std::setw(10) << n 
                  << std::setw(15) << std::fixed << std::setprecision(6) << final_grid->x[n]
                  << std::setw(15) << u[n]
                  << std::setw(15) << v[n]
                  << std::setw(20) << std::scientific << (u[n] - v[n]) << "\n";
    }
    std::cout << std::string(75, '-') << "\n";

    delete final_grid;
    return 0;
}