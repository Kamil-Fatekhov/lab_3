add_executable(p3_mixed_test_cls 
    main.cpp
    solver.cpp
)

target_include_directories(p3_mixed_test_cls PRIVATE
    ${CMAKE_SOURCE_DIR}/shared
)

target_link_libraries(p3_mixed_test_cls PRIVATE
    shared
    $<$<PLATFORM_ID:Linux,Darwin>:m>
)