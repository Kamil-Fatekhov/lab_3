#include <iostream>
#include <iomanip>
#include <vector>
#include "solver.h"

// Функция для формирования JSON вывода
void print_json(const BVP_Result& res) {
    std::cout << "{\n";
    std::cout << "  \"n\": " << res.grid.n << ",\n";
    std::cout << "  \"data\": [\n";
    for (size_t i = 0; i < res.u.size(); ++i) {
        std::cout << "    {\"x\": " << res.grid.x[i] << ", \"v\": " << res.u[i] << "}" 
                  << (i == res.u.size() - 1 ? "" : ",") << "\n";
    }
    std::cout << "  ]\n";
    std::cout << "}" << std::endl;
}

int main(int argc, char* argv[]) {
    int n;
    
    // Если n передано аргументом: ./main_imp.exe 1000
    if (argc > 1) {
        n = std::atoi(argv[1]);
    } else {
        std::cout << "Введите число разбиений n: ";
        std::cin >> n;
    }

    // Решаем задачу
    BVP_Result res = solve_mixed_main_imp(n);

    // Выводим данные в формате JSON для фронтенда
    print_json(res);

    return 0;
}