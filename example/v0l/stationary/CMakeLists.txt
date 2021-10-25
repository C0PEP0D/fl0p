cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(stationary)

# Build target
add_executable(stationary main.cpp)
target_include_directories(stationary PRIVATE .)

# Link modules
target_include_directories(stationary PRIVATE ../../../include)
target_include_directories(stationary PRIVATE ../../../../fl0w/include)
target_include_directories(stationary PRIVATE ../../../../p0l/include)
target_include_directories(stationary PRIVATE ../../../../v0l/include)
target_include_directories(stationary PRIVATE ../../../../m0sh/include)
# Link thirdparty
target_include_directories(stationary PRIVATE ../../../../../thirdparty/eigen)
