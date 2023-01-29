# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED "ON")

# Create project for serial jacobi
project(SERIAL_JACOBI LANGUAGES CXX)
add_executable(SERIAL_JACOBI src/jacobi_demo.cpp)

# Create project for parallel jacobi
project(PARALLEL_JACOBI LANGUAGES CXX)
add_executable(PARALLEL_JACOBI src/jacobi_demo.cpp)

# Create project for serial multigrid methods
project(SERIAL_MG_METHODS LANGUAGES CXX)
add_executable(SERIAL_MG_METHODS src/mg_demo.cpp)

# Create project for parallel multigrid methods
project(PARALLEL_MG_METHODS LANGUAGES CXX)
add_executable(PARALLEL_MG_METHODS src/mg_demo.cpp)

# Create project for serial multilevel multigrid methods
project(SERIAL_ML_MG_METHODS LANGUAGES CXX)
add_executable(SERIAL_ML_MG_METHODS src/multilevel_mg_demo.cpp)

# Create project for parallel multilevel multigrid methods
project(PARALLEL_ML_MG_METHODS LANGUAGES CXX)
add_executable(PARALLEL_ML_MG_METHODS src/multilevel_mg_demo.cpp)

# Locate OpenMP package
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(PARALLEL_MG_METHODS PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(PARALLEL_JACOBI PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(PARALLEL_ML_MG_METHODS PRIVATE OpenMP::OpenMP_CXX)
endif()

add_compile_options(-Wall -Wextra)









