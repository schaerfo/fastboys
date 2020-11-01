cmake_minimum_required(VERSION 3.16)
project(fastboys)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED On)

find_package(Boost REQUIRED COMPONENTS container)
find_package(Eigen3 REQUIRED)
find_package(fmt REQUIRED)
find_package(nlohmann_json REQUIRED)

add_executable(fastboys src/main.cpp
                        src/Molecule.hpp
                        src/Molecule.cpp)
target_link_libraries(fastboys PRIVATE Boost::container Eigen3::Eigen fmt::fmt nlohmann_json::nlohmann_json)
