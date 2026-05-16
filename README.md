# lab3 — Краевые задачи для ОДУ

## Структура проекта

```
lab3/
├── CMakeLists.txt          корневой — собирает всё
├── shared/                 общие утилиты (только заголовки)
│   ├── tridiag.h           метод прогонки
│   └── grid.h              сетка + интегралы для метода баланса
└── core/
    ├── p1_first_test/      чел. 1
    ├── p2_first_main/      чел. 2
    ├── p3_mixed_test_cls/  чел. 3
    ├── p4_mixed_main_imp/  чел. 4
    └── p5_mixed_test_imp/  чел. 5
```

---

## Как написать свой CMakeLists.txt

Скопируй этот шаблон в свою папку `core/pN_.../CMakeLists.txt`:

```cmake
add_executable(p2_first_main   # <- замени на своё имя папки
    main.cpp
    solver.cpp
)

target_link_libraries(p2_first_main PRIVATE   # <- то же имя
    shared
    $<$<PLATFORM_ID:Linux,Darwin>:m>
)
```

Больше ничего добавлять не нужно. `shared` автоматически пробрасывает
пути к `tridiag.h` и `grid.h` — дополнительный `target_include_directories`
не нужен.

После создания файла добавь свою строку в `core/CMakeLists.txt`:

```cmake
add_subdirectory(p2_first_main)   # <- своя папка
```

---

## Как подключить shared в коде

```cpp
#include "../../shared/tridiag.h"
#include "../../shared/grid.h"
```

---

## Сборка (WSL / Linux / macOS)

```bash
# Из корня проекта — один раз
cmake -B build -DCMAKE_BUILD_TYPE=Release

# Пересборка после изменений
cmake --build build

# Собрать только свою цель
cmake --build build --target p2_first_main

# Запуск
./build/core/p2_first_main/p2_first_main
```

Повторять `cmake -B build` нужно только если менял CMakeLists.txt.
После правки `.cpp`/`.h` достаточно `cmake --build build`.

## Сборка (Windows — PowerShell)

```powershell
# Добавить cmake в PATH если не находит
$env:PATH += ";C:\Program Files\CMake\bin"

cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

---

## Минимальные заглушки для старта

Чтобы проект собирался пока задача не реализована:

**solver.cpp**
```cpp
// TODO
```

**main.cpp**
```cpp
int main() { return 0; }
```

---

## Что есть в shared

### tridiag.h — namespace `tridiag`

| Функция | Что делает |
|---|---|
| `solve(a, b, c, d)` | Прогонка, система `a*x[i-1] + b*x[i] + c*x[i+1] = d[i]` |
| `solve_balance(A,B,C,D,mu1,mu2)` | Прогонка для схемы `-A*v[i-1] + B*v[i] - C*v[i+1] = D[i]` с ГУ первого рода |

### grid.h — namespace `grid`

| Что | Описание |
|---|---|
| `UniformGrid(n)` | Сетка: `x[i]=i*h`, `h=1/n`. Методы: `half(i)`, `left_half(i)`, `right_half(i)` |
| `integrate_simpson(f, a, b)` | Интеграл функции на отрезке, метод Симпсона |
| `k_half(k, xi, xi1)` | `(1/h)∫k dx` на полуинтервале — коэффициент A/C схемы |
| `q_bar(q, xl, xr)` | `(1/h)∫q dx` по ячейке — слагаемое в B |
| `f_bar(f, xl, xr)` | `(1/h)∫f dx` по ячейке — правая часть D |
| `k_half_jump(k1,k2,...)` | То же, но ячейка содержит точку разрыва ξ |
| `q_bar_jump(...)` | Аналогично для q |
| `f_bar_jump(...)` | Аналогично для f |
| `jump_node(g, xi)` | Индекс узла m, чья ячейка содержит ξ |
| `max_norm(u, v)` | `max|u[i]-v[i]|` → возвращает `{значение, индекс}` |
| `max_norm_coarse(v1, v2)` | `max|v1[i]-v2[2i]|` для сеток n и 2n → ε₂ |