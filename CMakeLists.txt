# Set minimum required CMake version
cmake_minimum_required(VERSION 3.16.3)

# Set project name
project(thames)

# Set options
option(THAMES_USE_SMARTUQ "Use SMART-UQ" ON)

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# Set C++ compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -g -fprofile-arcs -ftest-coverage")

# Add source code directory
add_subdirectory(src)

# Add SMART-UQ submodule
if(THAMES_USE_SMARTUQ)
    # Force SMART-UQ Release Build
    set(CMAKE_BUILD_TYPE "Release")
    # Force SMART-UQ to not build main file
    set(BUILD_MAIN OFF CACHE INTERNAL "" FORCE)
    # Add SMART-UQ directory
    add_subdirectory(external/smart-uq external/smart-uq)
endif(THAMES_USE_SMARTUQ)

# Build main executable
add_executable(main src/main.cpp)
# Link thames library
target_link_libraries(main thames)
# Link SMART-UQ static library
if(THAMES_USE_SMARTUQ)
    target_link_libraries(main smart-uq)
endif(THAMES_USE_SMARTUQ)
