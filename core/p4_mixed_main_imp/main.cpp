#include <iostream>
#include <iomanip>
#include <vector>
#include "solver.h"

void print_json(const BVP_Result& res) {
    std::cout << std::setprecision(10);
    std::cout << "{\n";
    std::cout << "  \"n\": "         << res.grid.n          << ",\n";
    std::cout << "  \"h\": "         << res.grid.h          << ",\n";
    std::cout << "  \"x_min\": "     << res.grid.x.front()  << ",\n";
    std::cout << "  \"x_max\": "     << res.grid.x.back()   << ",\n";
    std::cout << "  \"xi_jump\": "   << res.xi_jump         << ",\n";
    std::cout << "  \"gamma1\": "    << res.gamma1          << ",\n";
    std::cout << "  \"theta1\": "    << res.theta1          << ",\n";
    std::cout << "  \"gamma2\": "    << res.gamma2          << ",\n";
    std::cout << "  \"theta2\": "    << res.theta2          << ",\n";
    std::cout << "  \"bc_left\": \""  << res.bc_left  << "\",\n";
    std::cout << "  \"bc_right\": \"" << res.bc_right << "\",\n";
    std::cout << "  \"data\": [\n";
    for (size_t i = 0; i < res.u.size(); ++i) {
        std::cout << "    {\"x\": " << res.grid.x[i]
                  << ", \"u\": "    << res.u[i] << "}"
                  << (i == res.u.size() - 1 ? "" : ",") << "\n";
    }
    std::cout << "  ]\n";
    std::cout << "}" << std::endl;
}

int main(int argc, char* argv[]) {
    int n;
    if (argc > 1) {
        n = std::atoi(argv[1]);
    } else {
        std::cout << "Введите число разбиений n: ";
        std::cin >> n;
    }

    BVP_Result res = solve_mixed_main_imp(n);
    print_json(res);
    return 0;
}