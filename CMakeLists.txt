project(LP_MP-horizon-tracking)
cmake_minimum_required(VERSION 2.8.12)

# C++11
add_compile_options(-std=c++17)

# compiler options
add_definitions(-march=native)

add_subdirectory(external/LP_MP-MRF)

add_library(LP_MP-horizon-tracking INTERFACE)
target_include_directories(LP_MP-horizon-tracking INTERFACE include)
target_link_libraries(LP_MP-horizon-tracking INTERFACE LP_MP-MRF)

add_subdirectory(src)
enable_testing()
add_subdirectory(test)
