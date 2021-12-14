#ifndef TYPES
#define TYPES

#include <functional>

#include <Eigen/Core>

namespace thames::types {

    // Define vector types
    typedef Eigen::Matrix<double, 3, 1> Vector3;
    typedef Eigen::Matrix<double, 6, 1> Vector6;

    // Define matrix types
    typedef Eigen::Matrix<double, 3, 3> Matrix33;

    // Function types
    typedef std::function<double (double, Vector3)> Potential;
    typedef std::function<double (double, Vector3, Vector3)> PotentialDerivative;
    typedef std::function<Vector3 (double, Vector3, Vector3)> Force;
    
    // Dimensional factors structure
    typedef struct {
        double time;
        double length;
        double velocity;
        double grav;
    } DimensionalFactors;
    
}

#endif