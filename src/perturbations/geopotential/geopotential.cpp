#include <iostream>

#include "geopotential.h"
#include "../../types.h"

using namespace thames::types;

namespace thames::perturbations::geopotential{

    J2::J2(double mu, double J2, double radius){
        m_mu = mu;
        m_J2 = J2;
        m_radius = radius;
    }

    Vector3 J2::acceleration_total(double t, Vector3 R, Vector3 V){
        // Extract position components
        double x = R[0], y = R[1], z = R[2];

        // Calculate range
        double r = R.norm();

        // Precompute factors
        double J2_fac1 = -1.5*m_mu*m_J2*pow(m_radius, 2.0)/pow(r, 5.0);
        double J2_fac2 = 5.0*pow(z, 2.0)/pow(r, 2.0);

        // Declare and calculate perturbing acceleration vector
        Vector3 A;
        A << J2_fac1*x*(1.0 - J2_fac2),
             J2_fac1*y*(1.0 - J2_fac2),
             J2_fac1*z*(3.0 - J2_fac2);

        // Return perturbing acceleration vector
        return A;
    }

    double J2::potential(double t, Vector3 R){
        // Extract position components
        double x = R[0], y = R[1], z = R[2];

        // Calculate range
        double r = R.norm();

        // Calculate cosine of latitude
        double cphi = z/r;

        // Calculate perturbing potential
        double U = 0.5*m_J2*m_mu/pow(r, 3.0)*pow(m_radius, 2.0)*(3.0*pow(cphi, 2.0) - 1.0);

        // Return perturbing potential
        return U;
    }

}