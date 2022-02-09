if(THAMES_BUILD_STATIC)
    # Build thames static library
    add_library(${PROJECT_NAME} STATIC
                # Conversions
                conversions/cartesian.cpp
                conversions/dimensional.cpp
                conversions/geqoe.cpp
                conversions/keplerian.cpp
                # Perturbations
                perturbations/baseperturbation.cpp
                perturbations/geopotential/J2.cpp
                # Propagators
                propagators/cowell.cpp
                propagators/geqoe.cpp
                # Util
                util/angles.cpp
                util/root.cpp
                util/polynomials.cpp
                util/optimise.cpp
                util/sampling.cpp
                # Vector
                vector/arithmeticoverloads.cpp
                vector/geometry.cpp
    )
endif(THAMES_BUILD_STATIC)

if(THAMES_BUILD_MAIN)
    # Build main executable
    add_executable(main main.cpp)
    # Link thames library
    target_link_libraries(main thames)
    # Link SMART-UQ static library and FFTW libraries
    if(THAMES_USE_SMARTUQ)
        target_link_libraries(main smart-uq)
        target_link_libraries(main ${FFTW_LIBRARIES})
    endif(THAMES_USE_SMARTUQ)
endif(THAMES_BUILD_MAIN)
