cmake_minimum_required(VERSION 3.16.3)

project(lab2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(external/src/eigen)

set(EXECUTABLE_NAME main)

add_executable(${EXECUTABLE_NAME})

target_sources(${EXECUTABLE_NAME}
  PRIVATE
    src/main.cpp

    external/src/contract/src/contract.cpp
)

target_include_directories(${EXECUTABLE_NAME}
  PRIVATE
    include/public

    external/src/contract/include
)


target_link_libraries(${EXECUTABLE_NAME}
  PRIVATE
    Eigen3::Eigen
)
