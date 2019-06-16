cmake_minimum_required(VERSION 3.13)

# Helper functions for the unit tests
set(ODEX_TEST_DIR "${PROJECT_SOURCE_DIR}/test")
file(GLOB ODEX_TEST_HEADERS "${ODEX_TEST_DIR}/*.hpp")

# Test Executables
add_executable(Test_Threading "Test_Threading.cpp")
add_executable(Test_Partition "Test_Partition.cpp")
add_executable(Test_ExtrapolationStepper "Test_ExtrapolationStepper.cpp")

# Link required libraries for odex
odex_target_link_required_libraries(Test_Threading)
odex_target_link_required_libraries(Test_Partition)
odex_target_link_required_libraries(Test_ExtrapolationStepper)

# Eigen used in the ExtrapolationStepper test
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(Test_ExtrapolationStepper Eigen3::Eigen)

