#include "../types.h"

using namespace thames::types;

namespace thames::perturbations::geopotential{

    Vector3 J2_acceleration(double t, Vector3 R, Vector3 V, double mu, double J2, double radius){
        // TODO: documentation

        // Extract position components
        double x = R[0], y = R[1], z = R[2];

        // Calculate range
        double r = R.norm();

        // Precompute factors
        double J2_fac1 = -1.5*mu*J2*pow(radius, 2.0)/pow(r, 5.0);
        double J2_fac2 = 5.0*pow(z, 2.0)/pow(r, 2.0);

        // Declare and calculate perturbing acceleration vector
        Vector3 A;
        A << J2_fac1*x*(1.0 - J2_fac2),
             J2_fac1*y*(1.0 - J2_fac2),
             J2_fac1*z*(3.0 - J2_fac2);

        // Return perturbing acceleration vector
        return A;
    }

    double J2_potential(double t, Vector3 R, double mu, double J2, double radius){
        // TODO: documentation

        // Extract position components
        double x = R[0], y = R[1], z = R[2];

        // Calculate range
        double r = R.norm();

        // Calculate cosine of latitude
        double cphi = z/r;

        // Calculate perturbing potential
        double U = 0.5*J2*mu/pow(r, 3.0)*pow(radius, 2.0)*(3.0*pow(cphi, 2.0) - 1.0);

        // Return perturbing potential
        return U;
    }

    double J2_dpotential(double t, Vector3 R, Vector3 V){
        // TODO: documentation

        // Return zero change
        return 0.0;
    }

}