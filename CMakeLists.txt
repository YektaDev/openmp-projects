cmake_minimum_required(VERSION 3.30)
project(parallel)

cmake_policy(SET CMP0074 NEW) # Address the warning

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

add_executable(main main.cpp)

# Configure OpenMP for the target
if (OpenMP_FOUND)
    message(STATUS "OpenMP found")
    target_link_libraries(main OpenMP::OpenMP_CXX)
endif ()
