#ifndef THAMES_CONVERSIONS_DIMENSIONAL
#define THAMES_CONVERSIONS_DIMENSIONAL

#include <array>
#include <vector>

#include "../types.h"

using namespace thames::types;

namespace thames::conversions::dimensional{

    /**
     * @brief Non-dimensionalise Cartesian state.
     *
     * @tparam T Numeric type.
     * @param[in,out] t Current physical time.
     * @param[in,out] RV Current Cartesian state vector (position and velocity).
     * @param[in,out] mu Gravitational parameter.
     * @param[out] factors Structure containing the factors for non-dimensionalisation.
     */
    template<class T>
    void cartesian_nondimensionalise(T& t, std::array<T, 6>& RV, T& mu, DimensionalFactors& factors);

    /**
     * @brief Non-dimensionalise Cartesian state.
     *
     * @tparam T Numeric type.
     * @param[in,out] t Current physical time.
     * @param[in,out] RV Current Cartesian state vector (position and velocity).
     * @param[in,out] mu Gravitational parameter.
     * @param[out] factors Structure containing the factors for non-dimensionalisation.
     */
    template<class T>
    void cartesian_nondimensionalise(T& t, std::vector<T>& RV, T& mu, DimensionalFactors& factors);

    /**
     * @brief Dimensionalise Cartesian state.
     *
     * @tparam T Numeric type.
     * @param[in,out] t Current physical time.
     * @param[in,out] RV Current Cartesian state vector (position and velocity).
     * @param[in,out] mu Gravitational parameter.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     */ 
    template<class T>
    void cartesian_dimensionalise(T& t, std::array<T, 6>& RV, T& mu, const DimensionalFactors& factors);

    /**
     * @brief Dimensionalise Cartesian state.
     *
     * @tparam T Numeric type.
     * @param[in,out] t Current physical time.
     * @param[in,out] RV Current Cartesian state vector (position and velocity).
     * @param[in,out] mu Gravitational parameter.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     */ 
    template<class T>
    void cartesian_dimensionalise(T& t, std::vector<T>& RV, T& mu, const DimensionalFactors& factors);

}

#endif