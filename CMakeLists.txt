# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name and language
project(LightNucleiEmission CXX)

# Require C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set aggressive optimization flags for Release builds
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wextra")

# Find dependencies
find_package(Boost REQUIRED)
find_package(OpenMP REQUIRED)

# Define the executable and its source files
add_executable(lightnuclei 
    src/main.cpp 
    src/lightnuclei.cpp
)

# Include libraries
target_include_directories(lightnuclei PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(lightnuclei PRIVATE OpenMP::OpenMP_CXX)
