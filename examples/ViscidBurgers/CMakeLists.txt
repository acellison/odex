cmake_minimum_required(VERSION 3.13)

set(EXAMPLE_EXE Example_ViscidBurgers)
add_executable(${EXAMPLE_EXE} ViscidBurgers.cpp)

odex_example_link_libraries(${EXAMPLE_EXE})

