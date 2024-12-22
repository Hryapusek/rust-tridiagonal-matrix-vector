cmake_minimum_required(VERSION 3.16.3)

project(lab2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(external/src/eigen)

add_library(${PROJECT_NAME} STATIC
  external/src/contract/src/contract.cpp
    
  src/default_impl/euler_explicit_method.cc
  src/default_impl/euler_implicit_method.cc
  src/default_impl/main_matrix_calculator.cc
  src/default_impl/rkf_method.cc
    
  src/interval_splitter.cc
)

target_include_directories(${PROJECT_NAME}
  PUBLIC
    include/public

    external/src/contract/include
)

target_link_libraries(${PROJECT_NAME}
  PUBLIC
    Eigen3::Eigen
)

add_executable(${PROJECT_NAME}-main src/main.cc)
target_link_libraries(${PROJECT_NAME}-main ${PROJECT_NAME})
