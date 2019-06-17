cmake_minimum_required(VERSION 3.13)

set(EXAMPLE_EXE Example_KdV)
add_executable(${EXAMPLE_EXE} KdV.cpp)

odex_example_link_libraries(${EXAMPLE_EXE})

