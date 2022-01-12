#ifndef THAMES_CONVERSIONS_GEQOE
#define THAMES_CONVERSIONS_GEQOE

#include "../perturbations/baseperturbation.h"
#include "../types.h"

using namespace thames::types;
using namespace thames::perturbations::baseperturbation;

namespace thames::conversions::geqoe{

    /**
     * @brief Convert from Cartesian state to Generalised Equinoctial Orbital Elements (GEqOE).
     * 
     * @param[in] t Current physical time.
     * @param[in] RV Cartesian state.
     * @param[in] mu Gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @return Vector6 GEqOE state.
     */
    Vector6 cartesian_to_geqoe(const double &t, const Vector6 &RV, const double &mu, BasePerturbation<double, Vector3> &perturbation);

    /**
     * @brief Convert from Generalised Equinoctial Orbital Elements (GEqOE) to Cartesian state.
     * 
     * @param[in] t Current physical time.
     * @param[in] geqoe GEqOE state.
     * @param[in] mu Gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @return Vector6 Cartesian state.
     */
    Vector6 geqoe_to_cartesian(const double &t, const Vector6 &geqoe, const double &mu, BasePerturbation<double, Vector3> &perturbation);

}

#endif