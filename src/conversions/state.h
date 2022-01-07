#ifndef THAMES_CONVERSIONS_STATE
#define THAMES_CONVERSIONS_STATE

#include <Eigen/Core>

#include "../types.h"

using namespace thames::types;

namespace thames::conversions::state{

    /**
     * @brief Convert from Cartesian state to traditional Keplerian elements.
     * 
     * @param[in] RV Cartesian state.
     * @param[in] mu Gravitational parameter.
     * @return Vector6 Keplerian elements state.
     */
    Vector6 cartesian_to_keplerian(const Vector6 &RV, const double &mu);

    /**
     * @brief Convert from traditional Keplerian elements to Cartesian state.
     * 
     * @param[in] keplerian Keplerian elements state.
     * @param[in] mu Gravitational parameter.
     * @return Vector6 Cartesian state.
     */
    Vector6 keplerian_to_cartesian(const Vector6 &keplerian, const double &mu);

    /**
     * @brief Convert from Cartesian state to Generalised Equinoctial Orbital Elements (GEqOE).
     * 
     * @param[in] t Current physical time.
     * @param[in] RV Cartesian state.
     * @param[in] mu Gravitational parameter.
     * @param[in] U Perturbing potential function.
     * @return Vector6 GEqOE state.
     */
    Vector6 cartesian_to_geqoe(const double &t, const Vector6 &RV, const double &mu, const Potential &U);

    /**
     * @brief Convert from Generalised Equinoctial Orbital Elements (GEqOE) to Cartesian state.
     * 
     * @param[in] t Current physical time.
     * @param[in] geqoe GEqOE state.
     * @param[in] mu Gravitational parameter.
     * @param[in] U Perturbing potential function.
     * @return Vector6 Cartesian state.
     */
    Vector6 geqoe_to_cartesian(const double &t, const Vector6 &geqoe, const double &mu, const Potential &U);

}

#endif