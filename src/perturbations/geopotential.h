#ifndef THAMES_PERTURBATIONS_GEOPOTENTIAL
#define THAMES_PERTURBATIONS_GEOPOTENTIAL

#include "../types.h"

using namespace thames::types;

namespace thames::perturbations::geopotential{

    Vector3 J2_acceleration(double t, Vector3 R, Vector3 V, double mu, double J2, double radius);

    double J2_potential(double t, Vector3 R, double mu, double J2, double radius);

    double J2_dpotential(double t, Vector3 R, Vector3 V);

}

#endif