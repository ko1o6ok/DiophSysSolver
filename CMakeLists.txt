cmake_minimum_required(VERSION 3.21)
project(FastHNF)
set(CMAKE_CXX_STANDARD 14)

include_directories(include)
add_executable(FastHNF src/main.cpp src/Matrix.cpp src/SimplicialComplex.cpp src/Simplex.cpp)
