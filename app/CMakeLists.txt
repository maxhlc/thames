if(THAMES_BUILD_APP)
    # Add executable for GEqOE
    add_executable(geqoe geqoe.cpp)
    # Link executables
    target_link_libraries(geqoe thames)
    # Link SMART-UQ static library and FFTW libraries
    if(THAMES_USE_SMARTUQ)
        target_link_libraries(geqoe smart-uq)
        target_link_libraries(geqoe ${FFTW_LIBRARIES})
    endif(THAMES_USE_SMARTUQ)
endif(THAMES_BUILD_APP)