#include <array>
#include <cmath>

#include "dimensional.h"
#include "keplerian.h"
#include "../types.h"

using namespace thames::types;

namespace thames::conversions::dimensional{

    template<class T>
    void cartesian_nondimensionalise(T& t, std::array<T, 6>& RV, T& mu, DimensionalFactors& factors){
        // Calculate length factors (semi-major axis)
        std::array<T, 6> keplerian = thames::conversions::keplerian::cartesian_to_keplerian<T>(RV, mu);
        factors.length = keplerian[0];

        // Calculate velocity factor
        factors.velocity = sqrt(mu/factors.length);

        // Calculate time factors
        factors.time = sqrt(pow(factors.length, 3.0)/mu);

        // Calculate gravitational parameter factors
        factors.grav = mu;

        // Calculate non-dimensional states
        t /= factors.time;
        mu /= factors.grav;
        RV[0] /= factors.length;
        RV[1] /= factors.length;
        RV[2] /= factors.length;
        RV[3] /= factors.velocity;
        RV[4] /= factors.velocity;
        RV[5] /= factors.velocity;
    }
    template void cartesian_nondimensionalise<double>(double&, std::array<double, 6>&, double&, DimensionalFactors&);

    template<class T>
    void cartesian_dimensionalise(T& t, std::array<T, 6>& RV, T& mu, const DimensionalFactors& factors){
        // Calculate dimensional states
        t *= factors.time;
        mu *= factors.grav;
        RV[0] *= factors.length;
        RV[1] *= factors.length;
        RV[2] *= factors.length;
        RV[3] *= factors.velocity;
        RV[4] *= factors.velocity;
        RV[5] *= factors.velocity;
    }
    template void cartesian_dimensionalise<double>(double&, std::array<double, 6>&, double&, const DimensionalFactors&);

}