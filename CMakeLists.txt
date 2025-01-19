cmake_minimum_required(VERSION 3.30)
project(parallel)

set(CMAKE_CXX_STANDARD 26)

find_package(OpenMP REQUIRED)
if(OpenMP_FOUND)
    message(STATUS "OpenMP found!")
    message(STATUS "OpenMP include directories: ${OpenMP_CXX_INCLUDE_DIRS}")
    message(STATUS "OpenMP libraries: ${OpenMP_CXX_LIBRARIES}")
    message(STATUS "OpenMP compiler flags: ${OpenMP_CXX_FLAGS}")
else()
    message(FATAL_ERROR "OpenMP not found!")
endif()

add_executable(parallel main.cpp)

target_link_libraries(main OpenMP::OpenMP_CXX)
