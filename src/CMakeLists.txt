# Build thames static library
add_library(${PROJECT_NAME} STATIC
            # Conversions
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
            util/optimise.cpp
            # Vector
            vector/arithmeticoverloads.cpp
            vector/geometry.cpp
)
