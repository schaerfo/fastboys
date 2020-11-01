cmake_minimum_required(VERSION 3.16)
project(fastboys)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
find_package(fmt REQUIRED)

add_executable(fastboys src/main.cpp
                        src/Molecule.hpp
                        src/Molecule.cpp)
target_link_libraries(fastboys PRIVATE Eigen3::Eigen fmt::fmt)
