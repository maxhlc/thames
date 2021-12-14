#ifndef CONVERSIONS_STATE
#define CONVERSIONS_STATE

#include <Eigen/Core>

#include "../types.h"

using namespace thames::types;

namespace thames::conversions::state{

    Vector6 cartesian_to_keplerian(const Vector6 &RV, const double &mu);

    Vector6 keplerian_to_cartesian(const Vector6 &keplerian, const double &mu);

    Vector6 cartesian_to_geqoe(const double &t, const Vector6 &RV, const double &mu, const std::function<double (double, Vector3)> &U);

    Vector6 geqoe_to_cartesian(const double &t, const Vector6 &geqoe, const double &mu, const std::function<double (double, Vector3)> &U);

}

#endif