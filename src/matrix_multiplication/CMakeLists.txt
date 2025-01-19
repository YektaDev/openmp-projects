cmake_minimum_required(VERSION 3.30)
project(Parallel)

find_package(OpenMP REQUIRED)

add_executable(main main.cpp)

target_compile_features(main PUBLIC cxx_std_23)
target_link_libraries(main OpenMP::OpenMP_CXX)
