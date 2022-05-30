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

# Bundle source files together
# TODO: cleaner implementation
set(THAMES_SOURCE
    # Conversions
    conversions/cartesian.cpp
    conversions/dimensional.cpp
    conversions/geqoe.cpp
    conversions/keplerian.cpp
    conversions/universal.cpp
    # Input/output
    io/json.cpp
    io/point.cpp
    io/polynomial.cpp
    # Perturbations
    perturbations/atmosphere/atmospheremodel.cpp
    perturbations/atmosphere/drag.cpp
    perturbations/geopotential/J2.cpp
    perturbations/baseperturbation.cpp
    perturbations/perturbationcombiner.cpp
    # Propagators
    propagators/basepropagator.cpp
    propagators/cowell.cpp
    propagators/geqoe.cpp
    # Util
    util/angles.cpp
    util/optimise.cpp
    util/polynomials.cpp
    util/root.cpp
    util/sampling.cpp
    # Vector
    vector/arithmeticoverloads.cpp
    vector/geometry.cpp
)

# Bundle header files together
# TODO: cleaner implementation
set(THAMES_HEADERS
    # Constants
    ../include/constants/constants.h
    ../include/constants/earth.h
    ../include/constants/statetypes.h
    # Conversions
    ../include/conversions/cartesian.h
    ../include/conversions/conversions.h
    ../include/conversions/dimensional.h
    ../include/conversions/geqoe.h
    ../include/conversions/keplerian.h
    ../include/conversions/universal.h
    # Input/output
    ../include/io/io.h
    ../include/io/json.h
    ../include/io/point.h
    ../include/io/polynomial.h
    # Perturbations
    ../include/perturbations/atmosphere/atmospheremodel.h
    ../include/perturbations/atmosphere/drag.h
    ../include/perturbations/atmosphere/ussa76.h
    ../include/perturbations/geopotential/J2.h
    ../include/perturbations/baseperturbation.h
    ../include/perturbations/perturbationcombiner.h
    ../include/perturbations/perturbations.h
    # Propagators
    ../include/propagators/basepropagator.h
    ../include/propagators/cowell.h
    ../include/propagators/geqoe.h
    ../include/propagators/options.h
    ../include/propagators/propagators.h
    # Settings
    ../include/settings/settings.h
    # Util
    ../include/util/angles.h
    ../include/util/optimise.h
    ../include/util/polynomials.h
    ../include/util/root.h
    ../include/util/sampling.h
    ../include/util/util.h
    # Vector
    ../include/vector/arithmeticoverloads.h
    ../include/vector/geometry.h
    ../include/vector/vector.h

    # THAMES
    ../include/thames.h
)

if(THAMES_BUILD_STATIC)
    # Build thames static library
    add_library(${PROJECT_NAME} STATIC ${THAMES_SOURCE} ${THAMES_HEADERS})
    # Link to SMART-UQ
    if(THAMES_USE_SMARTUQ)
        target_link_libraries(${PROJECT_NAME} PUBLIC smart-uq)
    endif(THAMES_USE_SMARTUQ)
    # Link to nlohmann_json
    target_link_libraries(${PROJECT_NAME} PUBLIC nlohmann_json::nlohmann_json)
endif(THAMES_BUILD_STATIC)

if(THAMES_BUILD_MAIN)
    # Build main executable
    add_executable(main main.cpp)
    # Link thames library
    target_link_libraries(main ${PROJECT_NAME})
endif(THAMES_BUILD_MAIN)
