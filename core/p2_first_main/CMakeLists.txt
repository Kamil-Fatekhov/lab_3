add_executable(p2_first_main 
    main.cpp
    solver.cpp
)

target_include_directories(p2_first_main PRIVATE
    ${CMAKE_SOURCE_DIR}/shared
)

target_link_libraries(p2_first_main PRIVATE
    shared
    $<$<PLATFORM_ID:Linux,Darwin>:m>
)