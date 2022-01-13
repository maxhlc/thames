#ifndef THAMES_TYPES
#define THAMES_TYPES

#include <functional>

#include <Eigen/Core>

namespace thames::types {

    // Define vector types
    /// Type for a column vector of doubles (3 x 1)
    typedef Eigen::Matrix<double, 3, 1> Vector3;
    /// Type for a column vector of doubles (6 x 1)
    typedef Eigen::Matrix<double, 6, 1> Vector6;
    
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