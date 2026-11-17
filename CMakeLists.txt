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
target_include_directories(lightnuclei PRIVATE ${Boost_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(lightnuclei PRIVATE OpenMP::OpenMP_CXX)

# Optional test of the local momentum anisotropy (yield conservation, rest-frame moments):  cmake -DBUILD_TESTS=ON ..
option(BUILD_TESTS "Build tests/test_conservation" OFF)
if(BUILD_TESTS)
    add_executable(test_conservation tests/test_conservation.cpp src/lightnuclei.cpp)
    target_include_directories(test_conservation PRIVATE ${Boost_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(test_conservation PRIVATE OpenMP::OpenMP_CXX)
endif()
