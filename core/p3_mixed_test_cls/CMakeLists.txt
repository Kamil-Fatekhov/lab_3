add_library(p3_mixed_test_cls_lib main.cpp solver.cpp)

target_include_directories(p3_mixed_test_cls_lib PRIVATE
    ${CMAKE_SOURCE_DIR}/shared
)

target_link_libraries(p3_mixed_test_cls_lib PRIVATE
    shared
    $<$<PLATFORM_ID:Linux,Darwin>:m>
)