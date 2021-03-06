cmake_minimum_required(VERSION 3.16)
project(fastboys LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)

find_package(Boost REQUIRED COMPONENTS container program_options)
find_package(Eigen3 REQUIRED)
find_package(fmt REQUIRED)
find_package(nlohmann_json REQUIRED)

if (NOT FASTBOYS_DISABLE_OPENMP)
    find_package(OpenMP)
endif ()

if(MSVC)
    add_compile_options(/W3)
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

add_library(fastboys_common OBJECT src/Molecule.hpp
                                   src/Molecule.cpp
                                   src/OneElectronIntegrals.hpp
                                   src/OneElectronIntegralDetail.hpp
                                   src/OneElectronIntegrals.cpp
                                   src/TwoElectronIntegral.hpp
                                   src/TwoElectronIntegral.cpp
                                   src/FockMatrix.hpp
                                   src/FockMatrix.cpp
                                   src/DensityMatrix.cpp
                                   src/DensityMatrix.hpp
                                   src/Helpers.hpp)

target_link_libraries(fastboys_common PUBLIC Boost::container Eigen3::Eigen fmt::fmt nlohmann_json::nlohmann_json)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_link_libraries(fastboys_common PUBLIC pthread)
endif()

if (OpenMP_FOUND)
    message(STATUS "Building with OpenMP")
    target_compile_definitions(fastboys_common PUBLIC ENABLE_OPENMP)
    target_link_libraries(fastboys_common PUBLIC OpenMP::OpenMP_CXX)
else()
    message(STATUS "OpenMP disabled or not found, limited multi-threading support")
endif()

add_executable(fastboys src/main.cpp)
target_link_libraries(fastboys PRIVATE Boost::program_options fastboys_common)

set(GTEST_PATH /usr/src/googletest/googletest)
file(GLOB TEST_SRC test/*.cpp)
find_package(GTest)
if (${GTEST_FOUND})
    add_executable(fastboys_test ${TEST_SRC})
    target_include_directories(fastboys_test PRIVATE src)
    target_link_libraries(fastboys_test fastboys_common GTest::Main)
    set(BUILD_TESTS On)
elseif(EXISTS ${GTEST_PATH})
    add_executable(fastboys_test ${TEST_SRC} ${GTEST_PATH}/src/gtest-all.cc ${GTEST_PATH}/src/gtest_main.cc)
    target_include_directories(fastboys_test SYSTEM PRIVATE ${GTEST_PATH}/include ${GTEST_PATH})
    target_include_directories(fastboys_test PRIVATE src)
    target_link_libraries(fastboys_test fastboys_common pthread)
    set(BUILD_TESTS On)
else()
    set(BUILD_TESTS Off)
endif()

if(BUILD_TESTS)
    enable_testing()
    gtest_discover_tests(fastboys_test)
    message(STATUS "Building tests")
else()
    message(STATUS "Googletest not found, not building tests")
endif()
