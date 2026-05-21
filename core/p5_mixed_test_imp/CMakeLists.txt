add_executable(p5_mixed_test_imp
    main.cpp
    solver.cpp
)

target_include_directories(p5_mixed_test_imp PRIVATE
    ${CMAKE_SOURCE_DIR}/shared
)

target_link_libraries(p5_mixed_test_imp PRIVATE   # <- то же имя
    shared
    $<$<PLATFORM_ID:Linux,Darwin>:m>
)