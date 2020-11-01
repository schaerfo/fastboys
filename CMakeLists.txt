cmake_minimum_required(VERSION 3.12)
project(fastboys)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)

find_package(Boost REQUIRED COMPONENTS container)
find_package(Eigen3 REQUIRED)
find_package(fmt REQUIRED)
find_package(nlohmann_json REQUIRED)

if(MSVC)
    add_compile_options(/W3)
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

add_executable(fastboys src/main.cpp
                        src/Molecule.hpp
                        src/Molecule.cpp)
target_link_libraries(fastboys PRIVATE Boost::container Eigen3::Eigen fmt::fmt nlohmann_json::nlohmann_json)
