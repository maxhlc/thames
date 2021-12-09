#include <array>

#include <Eigen/Core>

#include "../types.h"

using namespace thames::types;

namespace thames::conversions::state{

    Vector6 cartesian_to_keplerian(const Vector6 &RV, const double &mu);

    Vector6 keplerian_to_cartesian(const Vector6 &keplerian, const double &mu);

}