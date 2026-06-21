#include "httplib.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <filesystem>

// Подключаем заголовки. 
// ВАЖНО: Если возникнут ошибки "redefinition", 
// проверь, чтобы функции в .h файлах не дублировались.
#include "core/p1_first_test/solver.h"
#include "core/p2_first_main/solver.h"
#include "core/p3_mixed_test_cls/solver.h"
#include "core/p4_mixed_main_imp/solver.h"
#include "core/p5_mixed_test_imp/solver.h"

using namespace httplib;
namespace fs = std::filesystem;

// Вспомогательная функция для JSON
std::string to_json(const std::vector<double>& x, const std::vector<double>& v, const std::vector<double>& u, bool is_main_task) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(10);
    ss << "{ \"x\": [";
    for (size_t i = 0; i < x.size(); ++i) ss << x[i] << (i == x.size() - 1 ? "" : ",");
    ss << "], \"v\": [";
    for (size_t i = 0; i < v.size(); ++i) ss << v[i] << (i == v.size() - 1 ? "" : ",");
    ss << "], \"u\": [";
    for (size_t i = 0; i < u.size(); ++i) ss << u[i] << (i == u.size() - 1 ? "" : ",");
    ss << "], \"is_main\": " << (is_main_task ? "true" : "false") << " }";
    return ss.str();
}

int main() {
    Server svr;

    // Исправляем путь к index.html
    svr.Get("/", [](const Request&, Response& res) {
        std::string path = "frontend/index.html";
        if (!fs::exists(path)) {
            path = "../frontend/index.html"; // если запуск из build/Debug
        }
        
        std::ifstream ifs(path);
        if (ifs) {
            std::string content((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
            res.set_content(content, "text/html; charset=utf-8");
        } else {
            res.status = 404;
            res.set_content("Error: index.html not found. Check paths!", "text/plain");
        }
    });

    svr.Get("/solve", [](const Request& req, Response& res) {
        int id = req.has_param("id") ? std::stoi(req.get_param_value("id")) : 3;
        int n = req.has_param("n") ? std::stoi(req.get_param_value("n")) : 100;

        std::vector<double> x, v, u;
        bool is_main = false;

        if (id == 1) { // 1-й чел
            v = solve_p1_task(n);
            for(int i=0; i<=n; ++i) {
                double xi = (double)i/n;
                x.push_back(xi);
                u.push_back(get_analytical_p1(xi));
            }
        } 
        else if (id == 2) { // 2-й чел (Основная задача)
            is_main = true;
            v = solve_bvp(n);
            auto v2 = solve_bvp(2 * n); // Для основной задачи u - это решение на сетке 2n
            for(int i=0; i<=n; ++i) {
                x.push_back((double)i/n);
                u.push_back(v2[2*i]);
            }
        }
        else if (id == 3) { // 3-й чел (Ты)
            v = solve_test_mixed(n);
            for(int i=0; i<=n; ++i) {
                double xi = (double)i/n;
                x.push_back(xi);
                u.push_back(get_analytical_test(xi));
            }
        }
        else if (id == 4) { // 4-й чел (Основная задача)
            is_main = true;
            auto res_n = solve_mixed_main_imp(n);
            auto res_2n = solve_mixed_main_imp(2 * n);
            
            // Используем .u, как указано в solver.h 4-го человека
            v = res_n.u; 
            for(int i = 0; i <= n; ++i) {
                x.push_back((double)i / n);
                // Берем значения из более густой сетки (каждый второй узел)
                u.push_back(res_2n.u[2 * i]); 
            }
        }
        else if (id == 5) { // 5-й чел
            // ВНИМАНИЕ: Если здесь ошибка линковки, значит функции p3 и p5 называются одинаково.
            // Нужно в одном из файлов переименовать solve_test_mixed в solve_test_mixed_imp
            v = solve_test_mixed_p5(n); 
            for(int i=0; i<=n; ++i) {
                double xi = (double)i/n;
                x.push_back(xi);
                u.push_back(get_analytical_test_p5(xi));
            }
        }

        res.set_content(to_json(x, v, u, is_main), "application/json");
    });

    std::cout << "Server started at http://localhost:8080" << std::endl;
    svr.listen("0.0.0.0", 8080);
    return 0;
}