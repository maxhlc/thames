#include <array>

#include <eigen3/Eigen/Core>

#include "../types.h"

using namespace thames::types;

namespace thames::conversions::state{

    Vector6 cartesian_to_keplerian(const Vector3 &R, const Vector3 &V, const double &mu);

}