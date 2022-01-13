#include <cmath>
#include <iostream>

#include "J2.h"
#include "../../util/vector.h"
#include "../../types.h"

using namespace thames::types;

namespace thames::perturbations::geopotential{

    template <class real, class vector>
    J2<real, vector>::J2(real mu, real J2, real radius){
        m_mu = mu;
        m_J2 = J2;
        m_radius = radius;
    }

    template <class real, class vector>
    vector J2<real, vector>::acceleration_total(real t, vector R, vector V){
        // Extract position components
        real x = R[0], y = R[1], z = R[2];

        // Calculate range
        real r = thames::util::vector::norm3<real, vector>(R);

        // Precompute factors
        real J2_fac1 = -1.5*m_mu*m_J2*pow(m_radius, 2.0)/pow(r, 5.0);
        real J2_fac2 = 5.0*pow(z, 2.0)/pow(r, 2.0);

        // Declare and calculate perturbing acceleration vector
        vector A;
        A[0] =  J2_fac1*x*(1.0 - J2_fac2),
        A[1] =  J2_fac1*y*(1.0 - J2_fac2),
        A[2] =  J2_fac1*z*(3.0 - J2_fac2);

        // Return perturbing acceleration vector
        return A;
    }

    template <class real, class vector>
    real J2<real, vector>::potential(real t, vector R){
        // Extract position components
        real z = R[2];

        // Calculate range
        real r = thames::util::vector::norm3<real, vector>(R);

        // Calculate cosine of latitude
        real cphi = z/r;

        // Calculate perturbing potential
        real U = 0.5*m_J2*m_mu/pow(r, 3.0)*pow(m_radius, 2.0)*(3.0*pow(cphi, 2.0) - 1.0);

        // Return perturbing potential
        return U;
    }
    template class J2<double, Vector3>;

}