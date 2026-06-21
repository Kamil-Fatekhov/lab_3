cmake_minimum_required(VERSION 3.16)
project(lab3 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wextra>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

find_package(MPI)

add_subdirectory(shared)
add_subdirectory(core)

# ====================================================================
# СБОРКА ОБЩЕГО СЕРВЕРА ДЛЯ ФРОНТЕНДА
# ====================================================================

# Создаем главный исполняемый файл сервера в корне проекта
add_executable(lab3_server server_main.cpp)
if(MSVC)
    set_target_properties(lab3_server PROPERTIES LINK_FLAGS "/FORCE:MULTIPLE")
endif()

# Указываем серверу, где искать заголовочные файлы (shared)
target_include_directories(lab3_server PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/shared
)

# Линкуем библиотеки всех 5 участников к серверу
target_link_libraries(lab3_server PRIVATE
    shared
    p1_first_test_lib
    p2_first_main_lib
    p3_mixed_test_cls_lib
    p4_mixed_main_imp_lib
    p5_mixed_test_imp_lib
)

# Если у первого человека нашелся MPI, подключаем его и к серверу, чтобы не было ошибок линковки
if(MPI_CXX_FOUND)
    target_link_libraries(lab3_server PRIVATE MPI::MPI_CXX)
endif()