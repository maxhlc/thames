#ifndef THAMES_TYPES
#define THAMES_TYPES

#include <array>
#include <functional>

namespace thames::types {
    
    // Dimensional factors structure
    /// Structure to contain factors for (non)dimensionalisation of Cartesian states
    typedef struct {
        double time;
        double length;
        double velocity;
        double grav;
    } DimensionalFactors;

}

#endif