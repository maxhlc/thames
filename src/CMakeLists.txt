add_library(${PROJECT_NAME} STATIC
            conversions/dimensional.cpp
            conversions/geqoe.cpp
            conversions/keplerian.cpp
            perturbations/baseperturbation.cpp
            perturbations/geopotential/J2.cpp
            propagators/cowell.cpp
            propagators/geqoe.cpp
            util/angles.cpp
            util/root.cpp
            util/optimise.cpp
            vector/arithmeticoverloads.cpp
            vector/geometry.cpp
)

add_executable(main main.cpp)

target_link_libraries(main thames)