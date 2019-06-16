cmake_minimum_required(VERSION 3.13)

# Eigen is required for the examples
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Link to python for matplotlib-cpp
find_package(Python3 COMPONENTS Interpreter Development)

# Find and link against python
if (Python3_EXECUTABLE)
  # Find out numpy's include path
  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c
            "from __future__ import print_function\ntry: import numpy; print(numpy.get_include(), end='')\nexcept:pass\n"
            OUTPUT_VARIABLE __numpy_path)
  # And the numpy version
  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c
            "from __future__ import print_function\ntry: import numpy; print(numpy.__version__, end='')\nexcept:pass\n"
    OUTPUT_VARIABLE __numpy_version)
else()
  message(STATUS "Python executable not found.")
endif(Python3_EXECUTABLE)

set(PYTHON_NUMPY_INCLUDE_DIR ${__numpy_path})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NumPy REQUIRED_VARS PYTHON_NUMPY_INCLUDE_DIR
                                        VERSION_VAR __numpy_version)

function(odex_target_link_python target)
    target_include_directories(${target} SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/dependencies/matplotlib-cpp)
    target_include_directories(${target} SYSTEM PRIVATE ${Python3_INCLUDE_DIRS})
    target_include_directories(${target} PRIVATE "${PYTHON_NUMPY_INCLUDE_DIR}")
    target_link_libraries(${target} ${Python3_LIBRARIES})
endfunction(odex_target_link_python)

function(odex_example_link_libraries target)
    target_link_libraries(${target} Eigen3::Eigen)
    odex_target_link_required_libraries(${target})
    odex_target_link_python(${target})
endfunction(odex_example_link_libraries target)


add_subdirectory(VanDerPolOscillator)
add_subdirectory(LorenzAttractor)
add_subdirectory(ViscidBurgers)
add_subdirectory(Transport1D)
add_subdirectory(WaveEquation)
add_subdirectory(PluckedString)

