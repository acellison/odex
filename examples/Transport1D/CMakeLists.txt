cmake_minimum_required(VERSION 3.13)

set(EXAMPLE_EXE Example_Transport1D)
add_executable(${EXAMPLE_EXE} Transport1D.cpp)

odex_example_link_libraries(${EXAMPLE_EXE})

