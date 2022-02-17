# MIT License
#
# Copyright (c) 2021-2022 Max Hallgarten La Casta
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Set minimum required CMake version
cmake_minimum_required(VERSION 3.16.3)

# Set project name
project(thames)

# Set options
option(THAMES_BUILD_MAIN "Build THAMES main" ON)
option(THAMES_BUILD_STATIC "Build THAMES static library" ON)
option(THAMES_BUILD_APP "Build THAMES applications" ON)
option(THAMES_USE_SMARTUQ "Use SMART-UQ" ON)
option(THAMES_USE_FFTW "Use FFTW" OFF)

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# Download/update git submodules
# Solution from: https://cliutils.gitlab.io/modern-cmake/chapters/projects/submodule.html
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

# Add SMART-UQ submodule
if(THAMES_USE_SMARTUQ)
    # Force SMART-UQ Release Build
    set(CMAKE_BUILD_TYPE "Release")
    # Add SMART-UQ directory
    add_subdirectory(external/smart-uq external/smart-uq)
    # Add flag to include SMART-UQ
    add_definitions(-DTHAMES_USE_SMARTUQ)
    # Set FFTW to true
    set(THAMES_USE_FFTW ON)
endif(THAMES_USE_SMARTUQ)

# Find FFTW libraries
if(THAMES_USE_FFTW)
    find_library(FFTW_LIBRARY fftw3)
    find_library(FFTWF_LIBRARY fftw3f)
    find_library(FFTWL_LIBRARY fftw3l)
    set(FFTW_LIBRARIES "${FFTW_LIBRARY}" "${FFTWF_LIBRARY}" "${FFTWL_LIBRARY}")
endif(THAMES_USE_FFTW)

# Set C++ compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3")

# Add source code directory
add_subdirectory(src)

# Add applications code directory
add_subdirectory(app)
