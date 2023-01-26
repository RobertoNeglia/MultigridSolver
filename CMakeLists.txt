# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.1)

# Create project AMG
project(EXEC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED "ON")

add_compile_options(-Wall -Wextra)

add_executable(EXEC src/main.cpp)

