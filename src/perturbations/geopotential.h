#ifndef THAMES_PERTURBATIONS_GEOPOTENTIAL
#define THAMES_PERTURBATIONS_GEOPOTENTIAL

#include "../types.h"

using namespace thames::types;

namespace thames::perturbations::geopotential{

    /**
     * @brief Calculate perturbing acceleration resulting from the J2 term.
     * 
     * @param t Current physical time.
     * @param R Position vector.
     * @param V Velocity vector.
     * @param mu Central body gravitational parameter.
     * @param J2 Central body J2 term.
     * @param radius Central body radius.
     * @return Vector3 Perturbing acceleration.
     */
    Vector3 J2_acceleration(double t, Vector3 R, Vector3 V, double mu, double J2, double radius);

    /**
     * @brief Calculate perturbing potential resulting from the J2 term.
     * 
     * @param t Current physical time.
     * @param R Position vector.
     * @param mu Central body gravitational parameter.
     * @param J2 Central body J2 term.
     * @param radius Central body radius.
     * @return double Perturbing potential.
     */
    double J2_potential(double t, Vector3 R, double mu, double J2, double radius);

    /**
     * @brief Calculate time derivative of the perturbing potential resulting from the J2 term.
     * 
     * @param t Current physical time.
     * @param R Position vector.
     * @param V Velocity vector.
     * @return double Time derivative of the perturbing potential.
     */
    double J2_dpotential(double t, Vector3 R, Vector3 V);

}

#endif