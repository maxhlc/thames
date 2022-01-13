#ifndef THAMES_CONVERSIONS_DIMENSIONAL
#define THAMES_CONVERSIONS_DIMENSIONAL

#include "../types.h"

using namespace thames::types;

namespace thames::conversions::dimensional{

    /**
    * @brief Non-dimensionalise Cartesian state.
    *
    * @tparam real Type for real numbers.
    * @tparam vector3 Type for the state vector slices.
    * @tparam vector6 Type for the state vector.
    * @param[in,out] t Current physical time.
    * @param[in,out] RV Current Cartesian state vector (position and velocity).
    * @param[in,out] mu Gravitational parameter.
    * @param[out] factors Structure containing the factors for non-dimensionalisation.
    */
    template<class real, class vector3, class vector6>
    void cartesian_nondimensionalise(real &t, vector6 &RV, real &mu, DimensionalFactors &factors);

    /**
    * @brief Dimensionalise Cartesian state.
    *
    * @tparam real Type for real numbers.
    * @tparam vector3 Type for the state vector slices.
    * @tparam vector6 Type for the state vector.
    * @param[in,out] t Current physical time.
    * @param[in,out] RV Current Cartesian state vector (position and velocity).
    * @param[in,out] mu Gravitational parameter.
    * @param[in] factors Structure containing the factors for dimensionalisation.
    */ 
    template<class real, class vector3, class vector6>
    void cartesian_dimensionalise(real &t, vector6 &RV, real &mu, const DimensionalFactors &factors);

}

#endif