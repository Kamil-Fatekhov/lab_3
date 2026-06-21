# core/p2_first_main/CMakeLists.txt
add_library(p2_first_main_lib main.cpp solver.cpp)

target_include_directories(p2_first_main_lib PRIVATE
    ${CMAKE_SOURCE_DIR}/shared
)

target_link_libraries(p2_first_main_lib PRIVATE
    shared
    $<$<PLATFORM_ID:Linux,Darwin>:m>
)

# ВОТ ЭТУ СТРОЧКУ НАДО ДОБАВИТЬ:
target_compile_definitions(p2_first_main_lib PRIVATE _USE_MATH_DEFINES)