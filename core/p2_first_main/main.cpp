#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
 
#include "solver.h"
#include "../../shared/grid.h"

// Функция для вычисления порядка сходимости p = log(eps2_prev / eps2_curr) / log(2)
double calculate_order(double eps_prev, double eps_curr) {
    if (eps_prev <= 0 || eps_curr <= 0) return 0.0;
    return std::log2(eps_prev / eps_curr);
}

static void run_main_problem() {
    const double EPS = 0.5e-6;
    int n_final = 0;
    auto [v, eps2] = solve_with_accuracy(EPS, 10, n_final);
    
    // ... ваш код вывода справки ...
    std::cout << "\n=== ИТОГОВЫЙ АНАЛИЗ СХОДИМОСТИ (Основная задача) ===\n";
    std::cout << std::left << std::setw(10) << "n" 
              << std::setw(20) << "eps2" 
              << std::setw(15) << "Порядок p" << "\n";
    std::cout << std::string(45, '-') << "\n";

    // Здесь можно дописать цикл, который прогонит n от 10 до n_final 
    // и выведет таблицу с порядком сходимости:
    double prev_eps = 0.0;
    for (int n = 10; n <= n_final; n *= 2) {
        auto v_n = solve_bvp(n);
        auto v_2n = solve_bvp(2 * n);
        auto [err, idx] = compute_eps2(v_n, v_2n);
        
        double p = (prev_eps > 0) ? calculate_order(prev_eps, err) : 0.0;
        
        std::cout << std::left << std::setw(10) << n 
                  << std::scientific << std::setw(20) << err 
                  << std::fixed << std::setprecision(2) << std::setw(15) << p << "\n";
        
        prev_eps = err;
    }
}

int main() {
    run_main_problem();
    return 0;
}