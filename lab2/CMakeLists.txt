cmake_minimum_required(VERSION 3.16.3)

project(lab2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(external/src/eigen)

set(EXECUTABLE_NAME main)

add_executable(${EXECUTABLE_NAME})

target_sources(${EXECUTABLE_NAME}
  PUBLIC
    external/src/contract/src/contract.cpp
    
    src/default_impl/euler_explicit_method.cc
    src/default_impl/euler_implicit_method.cc
    src/default_impl/main_matrix_calculator.cc
    
    src/interval_splitter.cc
    src/main.cc
)

target_include_directories(${EXECUTABLE_NAME}
  PUBLIC
    include/public

    external/src/contract/include
)


target_link_libraries(${EXECUTABLE_NAME}
  PUBLIC
    Eigen3::Eigen
)
