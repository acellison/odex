cmake_minimum_required(VERSION 3.13)

# C++ Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything -Wno-deprecated-register -Wno-c++98-compat -Wno-double-promotion -Wno-padded")

# odex include directories
set(ODEX_INCLUDE "${PROJECT_SOURCE_DIR}/include")

# Set up source groups for the odex library so we can view them in the IDE.
# The call to add_custom_target creates a dummy target that simply lists the
# sources in hierarchy mirroring that of the odex include directory structure.
file(GLOB ODEX_HEADERS_BASE "${ODEX_INCLUDE}/odex/*.hpp")
file(GLOB ODEX_HEADERS_STEPPERS "${ODEX_INCLUDE}/odex/steppers/*.hpp")
file(GLOB ODEX_HEADERS_OBSERVERS "${ODEX_INCLUDE}/odex/observers/*.hpp")
file(GLOB ODEX_HEADERS_THREADING "${ODEX_INCLUDE}/odex/threading/*.hpp")
file(GLOB ODEX_HEADERS_DETAIL "${ODEX_INCLUDE}/odex/detail/*.hpp")
source_group(odex FILES ${ODEX_HEADERS_BASE})
source_group(odex\\steppers FILES ${ODEX_HEADERS_STEPPERS})
source_group(odex\\observers FILES ${ODEX_HEADERS_OBSERVERS})
source_group(odex\\threading FILES ${ODEX_HEADERS_THREADING})
source_group(odex\\detail FILES ${ODEX_HEADERS_DETAIL})
add_custom_target(odex_sources SOURCES
    ${ODEX_HEADERS_BASE}
    ${ODEX_HEADERS_STEPPERS}
    ${ODEX_HEADERS_OBSERVERS}
    ${ODEX_HEADERS_THREADING}
    ${ODEX_HEADERS_DETAIL}
  )

# Threads package is required for odex::extrapolation_stepper
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Clients can use this to link to odex's required libraries
function(odex_target_link_required_libraries target)
    target_include_directories(${target} PRIVATE ${ODEX_INCLUDE})
    target_link_libraries(${target} Threads::Threads)
endfunction(odex_target_link_required_libraries target)


