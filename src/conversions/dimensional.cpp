#include "dimensional.h"
#include "state.h"
#include "../types.h"

using namespace thames::types;

namespace thames::conversions::dimensional{

    void cartesian_nondimensionalise(double &t, Vector6 &RV, double &mu, DimensionalFactors &factors){
        // Extract position and velocity vectors
        Vector3 R, V;
        R = RV(Eigen::seq(0,2));
        V = RV(Eigen::seq(3,5));

        // Calculate length factors (semi-major axis)
        Vector6 keplerian = thames::conversions::state::cartesian_to_keplerian(RV, mu);
        factors.length = keplerian[0];

        // Calculate velocity factor
        factors.velocity = sqrt(mu/factors.length);

        // Calculate time factors
        factors.time = sqrt(pow(factors.length, 3.0)/mu);

        // Calculate gravitational parameter factors
        factors.grav = mu;

        // Calculate non-dimensional states
        t = t/factors.time;
        mu = mu/factors.grav;
        R = R/factors.length;
        V = V/factors.velocity;
        RV << R, V;
    }

    void cartesian_dimensionalise(double &t, Vector6 &RV, double &mu, const DimensionalFactors &factors){
        // Extract position and velocity vectors
        Vector3 R, V;
        R = RV(Eigen::seq(0,2));
        V = RV(Eigen::seq(3,5));

        // Calculate dimensional states
        t = t*factors.time;
        mu = mu*factors.grav;
        R = R*factors.length;
        V = V*factors.velocity;
        RV << R, V;
    }

}