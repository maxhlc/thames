#ifndef THAMES_CONVERSIONS_GEQOE
#define THAMES_CONVERSIONS_GEQOE

#include <array>
#include <vector>

#include "../perturbations/baseperturbation.h"

using namespace thames::perturbations::baseperturbation;

namespace thames::conversions::geqoe{

    ////////////
    // Arrays //
    ////////////

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
    std::array<T, 6> cartesian_to_geqoe(const T& t, const std::array<T, 6>& RV, const T& mu, BasePerturbation<T>* perturbation);

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
    std::array<T, 6> geqoe_to_cartesian(const T& t, const std::array<T, 6>& geqoe, const T& mu, BasePerturbation<T>* perturbation);

    /////////////
    // Vectors //
    /////////////

    /**
     * @brief Convert from Cartesian state to Generalised Equinoctial Orbital Elements (GEqOE).
     * 
     * @tparam T Numeric type.
     * @param[in] t Current physical time.
     * @param[in] RV Cartesian state.
     * @param[in] mu Gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @return std::vector<T> GEqOE state.
     */
    template<class T>
    std::vector<T> cartesian_to_geqoe(const T& t, const std::vector<T>& RV, const T& mu, BasePerturbation<T>* perturbation);

    /**
     * @brief Convert from Generalised Equinoctial Orbital Elements (GEqOE) to Cartesian state.
     * 
     * @tparam T Numeric type.
     * @param[in] t Current physical time.
     * @param[in] geqoe GEqOE state.
     * @param[in] mu Gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @return std::vector<T> Cartesian state.
     */
    template<class T>
    std::vector<T> geqoe_to_cartesian(const T& t, const std::vector<T>& geqoe, const T& mu, BasePerturbation<T>* perturbation);

}

#endif