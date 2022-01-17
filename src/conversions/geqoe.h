#ifndef THAMES_CONVERSIONS_GEQOE
#define THAMES_CONVERSIONS_GEQOE

#include <array>

#include "../perturbations/baseperturbation.h"

using namespace thames::perturbations::baseperturbation;

namespace thames::conversions::geqoe{

    /**
     * @brief Convert from Cartesian state to Generalised Equinoctial Orbital Elements (GEqOE).
     * 
     * @tparam T Numeric type.
     * @param[in] t Current physical time.
     * @param[in] RV Cartesian state.
     * @param[in] mu Gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @return std::array<T, 6> GEqOE state.
     */
    template<class T>
    std::array<T, 6> cartesian_to_geqoe(const T& t, const std::array<T, 6>& RV, const T& mu, const BasePerturbation<T>& perturbation);

    /**
     * @brief Convert from Generalised Equinoctial Orbital Elements (GEqOE) to Cartesian state.
     * 
     * @tparam T Numeric type.
     * @param[in] t Current physical time.
     * @param[in] geqoe GEqOE state.
     * @param[in] mu Gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @return std::array<T, 6> Cartesian state.
     */
    template<class T>
    std::array<T, 6> geqoe_to_cartesian(const T& t, const std::array<T, 6>& geqoe, const T& mu, const BasePerturbation<T>& perturbation);

}

#endif