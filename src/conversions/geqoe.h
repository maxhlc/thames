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
     * @tparam real Type for real numbers.
     * @tparam vector3 Type for the state vector slices.
     * @tparam vector6 Type for the state vector.
     * @param[in] t Current physical time.
     * @param[in] RV Cartesian state.
     * @param[in] mu Gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @return Vector6 GEqOE state.
     */
    template<class real, class vector3, class vector6>
    vector6 cartesian_to_geqoe(const real &t, const vector6 &RV, const real &mu, BasePerturbation<real, vector3> &perturbation);

    /**
     * @brief Convert from Generalised Equinoctial Orbital Elements (GEqOE) to Cartesian state.
     * 
     * @tparam real Type for real numbers.
     * @tparam vector3 Type for the state vector slices.
     * @tparam vector6 Type for the state vector.
     * @param[in] t Current physical time.
     * @param[in] geqoe GEqOE state.
     * @param[in] mu Gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @return Vector6 Cartesian state.
     */
    template<class real, class vector3, class vector6>
    vector6 geqoe_to_cartesian(const real &t, const vector6 &geqoe, const real &mu, BasePerturbation<real, vector3> &perturbation);

}

#endif