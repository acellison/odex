cmake_minimum_required(VERSION 3.13)

set(EXAMPLE_EXE Example_PluckedString)
add_executable(${EXAMPLE_EXE} PluckedString.cpp)

odex_example_link_libraries(${EXAMPLE_EXE})

