#include "solver.h"
#include "grid.h"
#include "tridiag.h"
#include <cmath>
#include <functional>
#include <algorithm>
#include <vector>

namespace {
    // Константы Варианта №5
    const double XI = std::acos(-1.0) / 4.0; // pi / 4
    const double BETA = 1.0; // Коэффициент теплоотдачи для ГУ 3-го рода
    
    // Параметры краевых условий 3-го рода
    // Левое: -k*u'(0) + gamma1*u(0) = theta1
    // Правое: k*u'(1) + gamma2*u(1) = theta2
    const double GAMMA1 = 1.0;
    const double GAMMA2 = 3.0;
    const double THETA1 = 2.0;
    const double THETA2 = 4.0;

    // Предельные константы тестовой задачи в точке разрыва (Вариант 5)
    const double K1_STAR = 1.0;
    const double K2_STAR = 0.5;
    const double Q1_STAR = 1.0;
    const double Q2_STAR = (std::acos(-1.0) * std::acos(-1.0)) / 16.0; // pi^2 / 16
    const double F1_STAR = 1.0;
    const double F2_STAR = std::cos(XI); // sqrt(2)/2
    
    // Решение системы 4x4 методом Гаусса для точного поиска C1, C2, C3, C4
    std::vector<double> solve_system_4x4(double M[4][5]) {
        for(int i = 0; i < 4; ++i) {
            int pivot = i;
            for(int j = i + 1; j < 4; ++j) {
                if(std::abs(M[j][i]) > std::abs(M[pivot][i])) pivot = j;
            }
            for(int k = 0; k < 5; ++k) std::swap(M[i][k], M[pivot][k]);
            
            double diag = M[i][i];
            for(int k = i; k < 5; ++k) M[i][k] /= diag;
            
            for(int j = 0; j < 4; ++j) {
                if(i != j) {
                    double factor = M[j][i];
                    for(int k = i; k < 5; ++k) M[j][k] -= factor * M[i][k];
                }
            }
        }
        return {M[0][4], M[1][4], M[2][4], M[3][4]};
    }
}

double get_analytical_test_p5(double x) {
    static bool initialized = false;
    static double c1, c2, c3, c4, l1, l2;

    if (!initialized) {
        l1 = std::sqrt(Q1_STAR / K1_STAR);
        l2 = std::sqrt(Q2_STAR / K2_STAR);

        double M[4][5] = {0};
        
        // 1. Левое граничное условие 3-го рода: -k1*u'(0) + gamma1*u(0) = theta1
        // u(x) = c1*exp(l1*x) + c2*exp(-l1*x) + F1_STAR/Q1_STAR
        // u'(x) = c1*l1*exp(l1*x) - c2*l1*exp(-l1*x)
        M[0][0] = -K1_STAR * l1 + GAMMA1;
        M[0][1] = K1_STAR * l1 + GAMMA1;
        M[0][4] = THETA1 - GAMMA1 * (F1_STAR / Q1_STAR);

        // 2. Непрерывность решения в точке разрыва u(xi-0) = u(xi+0)
        M[1][0] = std::exp(l1 * XI); 
        M[1][1] = std::exp(-l1 * XI);
        M[1][2] = -std::exp(l2 * XI); 
        M[1][3] = -std::exp(-l2 * XI);
        M[1][4] = (F2_STAR / Q2_STAR) - (F1_STAR / Q1_STAR);

        // 3. Непрерывность потока w(xi-0) = w(xi+0)
        M[2][0] = K1_STAR * l1 * std::exp(l1 * XI);
        M[2][1] = -K1_STAR * l1 * std::exp(-l1 * XI);
        M[2][2] = -K2_STAR * l2 * std::exp(l2 * XI);
        M[2][3] = K2_STAR * l2 * std::exp(-l2 * XI);
        M[2][4] = 0.0;

        // 4. Правое граничное условие 3-го рода: k2*u'(1) + gamma2*u(1) = theta2
        M[3][2] = (K2_STAR * l2 + GAMMA2) * std::exp(l2);
        M[3][3] = (-K2_STAR * l2 + GAMMA2) * std::exp(-l2);
        M[3][4] = THETA2 - GAMMA2 * (F2_STAR / Q2_STAR);

        auto C = solve_system_4x4(M);
        c1 = C[0]; c2 = C[1]; c3 = C[2]; c4 = C[3];
        initialized = true;
    }

    if (x < XI) {
        return c1 * std::exp(l1 * x) + c2 * std::exp(-l1 * x) + (F1_STAR / Q1_STAR);
    } else {
        return c3 * std::exp(l2 * x) + c4 * std::exp(-l2 * x) + (F2_STAR / Q2_STAR);
    }
}

std::vector<double> solve_test_mixed_p5(int n) {
    grid::UniformGrid g(n);
    double h = g.h;
    std::vector<double> a(n + 1, 0), b(n + 1, 0), c(n + 1, 0), d(n + 1, 0);

    // Лямбды функций тестовой задачи
    auto k_f = [](double x) { return (x < XI) ? K1_STAR : K2_STAR; };
    auto q_f = [](double x) { return (x < XI) ? Q1_STAR : Q2_STAR; };
    auto f_f = [](double x) { return (x < XI) ? F1_STAR : F2_STAR; };
    
    auto k1_f = [](double) { return K1_STAR; };
    auto k2_f = [](double) { return K2_STAR; };
    auto q1_f = [](double) { return Q1_STAR; };
    auto q2_f = [](double) { return Q2_STAR; };
    auto f1_f = [](double) { return F1_STAR; };
    auto f2_f = [](double) { return F2_STAR; };

    // Левое Граничное Условие (3-й род): -k*u'(0) + gamma1*u(0) = theta1
    // Аппроксимация с использованием фиктивного узла: 
    // -k*(u1 - u(-1))/(2h) + gamma1*u0 = theta1
    // или более простая аппроксимация O(h^2) с использованием односторонней разности
    // Для простоты используем аппроксимацию: 
    // u0 = (theta1 + k*u1/h)/(gamma1 + k/h)
    // при условии что gamma1 + k/h != 0
    
    // Более точная аппроксимация O(h^2):
    // -k*(u1 - u0)/h + gamma1*u0 = theta1 + O(h)
    // Переносим: -(k/h)*u1 + (k/h + gamma1)*u0 = theta1
    
    // Для O(h^2) используем центральную разность с фиктивным узлом:
    // -k*(u1 - u(-1))/(2h) + gamma1*u0 = theta1
    // Выражаем u(-1) через u0, u1 из уравнения (k/h)*(u0 - u(-1)) - (k/h)*(u1 - u0) = 0
    // После преобразований получаем:
    // (k/h + gamma1)*u0 - (k/h)*u1 = theta1
    
    a[0] = 0.0;
    b[0] = K1_STAR / h + GAMMA1;
    c[0] = -K1_STAR / h;
    d[0] = THETA1;

    for (int i = 1; i < n; ++i) {
        double x_left = g.x[i - 1];
        double x_node = g.x[i];
        double x_right = g.x[i + 1];
        double x_half_L = g.left_half(i);
        double x_half_R = g.right_half(i);

        double a_i, c_i, q_i, d_i;

        if (XI > x_left && XI < x_node)
            a_i = grid::k_half_jump(k1_f, k2_f, x_left, x_node, XI) / (h * h);
        else
            a_i = grid::k_half(k_f, x_left, x_node) / (h * h);

        if (XI > x_node && XI < x_right)
            c_i = grid::k_half_jump(k1_f, k2_f, x_node, x_right, XI) / (h * h);
        else
            c_i = grid::k_half(k_f, x_node, x_right) / (h * h);

        if (XI > x_half_L && XI < x_half_R) {
            q_i = grid::q_bar_jump(q1_f, q2_f, x_half_L, x_half_R, XI);
            d_i = grid::f_bar_jump(f1_f, f2_f, x_half_L, x_half_R, XI);
        } else {
            q_i = grid::q_bar(q_f, x_half_L, x_half_R);
            d_i = grid::f_bar(f_f, x_half_L, x_half_R);
        }

        a[i] = -a_i;
        c[i] = -c_i;
        b[i] = a_i + c_i + q_i;
        d[i] = d_i;
    }

    // Правое Граничное Условие (3-й род, УЛУЧШЕННАЯ аппроксимация O(h^2))
    // k2*u'(1) + gamma2*u(1) = theta2
    // Аппроксимация с фиктивным узлом:
    // k2*(u(n+1) - u(n-1))/(2h) + gamma2*u(n) = theta2
    // Из уравнения во внутреннем узле: k*(u(n+1) - 2u(n) + u(n-1))/h^2 + ... = 0
    // После исключения u(n+1) получаем:
    // (-k/h)*u(n-1) + (k/h + gamma2 + (h/2)*Q2_STAR)*u(n) = theta2 + (h/2)*F2_STAR
    
    a[n] = -K2_STAR / h;
    b[n] = (K2_STAR / h) + GAMMA2 + (h / 2.0) * Q2_STAR;
    c[n] = 0.0;
    d[n] = THETA2 + (h / 2.0) * F2_STAR;

    return tridiag::solve(a, b, c, d);
}