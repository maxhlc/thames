#include <cmath>

#include "dimensional.h"
#include "keplerian.h"
#include "../types.h"

using namespace thames::types;

namespace thames::conversions::dimensional{

    template<class real, class vector3, class vector6>
    void cartesian_nondimensionalise(real &t, vector6 &RV, real &mu, DimensionalFactors &factors){
        // Calculate length factors (semi-major axis)
        vector6 keplerian = thames::conversions::keplerian::cartesian_to_keplerian<real, vector3, vector6>(RV, mu);
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
    template void cartesian_nondimensionalise<double, Vector3, Vector6>(double&, Vector6&, double&, DimensionalFactors&);

    template<class real, class vector3, class vector6>
    void cartesian_dimensionalise(real &t, vector6 &RV, real &mu, const DimensionalFactors &factors){
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
    template void cartesian_dimensionalise<double, Vector3, Vector6>(double&, Vector6&, double&, const DimensionalFactors&);

}