cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(unstationary)

# Build target
add_executable(unstationary main.cpp)
target_include_directories(unstationary PRIVATE .)

# Link modules
target_include_directories(unstationary PRIVATE ../../../include)
target_include_directories(unstationary PRIVATE ../../../../fl0w/include)
target_include_directories(unstationary PRIVATE ../../../../p0l/include)
target_include_directories(unstationary PRIVATE ../../../../v0l/include)
target_include_directories(unstationary PRIVATE ../../../../m0sh/include)
# Link thirdparty
target_include_directories(unstationary PRIVATE ../../../../../thirdparty/eigen)
