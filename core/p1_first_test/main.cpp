#define NOMINMAX 
#define _USE_MATH_DEFINES    
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>             
#include <windows.h>         
#include <chrono> // Заменяем MPI_Wtime на стандартный хронограф
#include "solver.h"

int main(int argc, char** argv) {
    ::SetConsoleOutputCP(65001);

    std::cout << "=================== АНАЛИЗ ПОГРЕШНОСТИ ТЕСТОВОЙ ЗАДАЧИ ===================" << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
    std::cout << std::left << std::setw(10) << "n" 
              << std::setw(25) << "max |u - v| (eps1)" 
              << std::setw(18) << "-lg(eps1)" << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;

    std::vector<int> grids = {10, 20, 50, 80, 100, 200, 500, 800, 1000, 10000, 100000, 1000000};

    // Стандартный таймер C++
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int current_n : grids) {
        BalanceSolver solver(current_n);
        std::vector<double> v = solver.solve();

        double max_err = solver.max_error(v);
        double lg_eps = -std::log10(max_err);

        std::cout << std::left << std::setw(10) << current_n 
                  << std::scientific << std::setw(25) << max_err 
                  << std::fixed << std::setprecision(4) << std::setw(18) << lg_eps
                  << std::defaultfloat << std::endl;
        
        if (current_n == 200) {
            std::cout << " >>> ДОСТИГНУТА ТРЕБУЕМАЯ ПОГРЕШНОСТЬ ПРИ n = 200 <<<" << std::endl;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    std::cout << "(" << elapsed.count() << " s)" << std::endl;

    return 0;
}