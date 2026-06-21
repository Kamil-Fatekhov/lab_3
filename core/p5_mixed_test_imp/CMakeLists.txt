add_library(p5_mixed_test_imp_lib main.cpp solver.cpp)
target_include_directories(p5_mixed_test_imp_lib PRIVATE ${CMAKE_SOURCE_DIR}/shared)
target_link_libraries(p5_mixed_test_imp_lib PRIVATE shared $<$<PLATFORM_ID:Linux,Darwin>:m>)