# core/p1_first_test/CMakeLists.txt
cmake_minimum_required(VERSION 3.10)

# Создаем интерфейсную (header-only) библиотеку, так как solver.cpp не существует
add_library(p1_first_test_lib INTERFACE)

# Говорим CMake, где серверу искать файлы math_utils.h и solver.h
target_include_directories(p1_first_test_lib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# Добавляем дефайны, которые автоматически применятся ко всем, кто подключит эту либу
target_compile_definitions(p1_first_test_lib INTERFACE _USE_MATH_DEFINES)