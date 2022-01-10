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

    // Define matrix types
    /// Type for a matrix of doubles (3 x 3)
    typedef Eigen::Matrix<double, 3, 3> Matrix33;

    // Function types
    /// Type for a function for potentials
    typedef std::function<double (double, Vector3)> PotentialFunc;
    /// Type for a function for time derivatives of potentials
    typedef std::function<double (double, Vector3, Vector3)> PotentialDerivativeFunc;
    /// Type for a function for accelerations
    typedef std::function<Vector3 (double, Vector3, Vector3)> AccelerationFunc;
    
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