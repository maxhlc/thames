#ifndef THAMES_CONVERSIONS_DIMENSIONAL
#define THAMES_CONVERSIONS_DIMENSIONAL

#include "../types.h"

using namespace thames::types;

namespace thames::conversions::dimensional{

    /**
    * @brief Non-dimensionalise Cartesian state.
    *
    * @param[in,out] t Current physical time.
    * @param[in,out] RV Current Cartesian state vector (position and velocity).
    * @param[in,out] mu Gravitational parameter.
    * @param[out] factors Structure containing the factors for non-dimensionalisation.
    */ 
    void cartesian_nondimensionalise(double &t, Vector6 &RV, double &mu, DimensionalFactors &factors);

    /**
    * @brief Dimensionalise Cartesian state.
    *
    * @param[in,out] t Current physical time.
    * @param[in,out] RV Current Cartesian state vector (position and velocity).
    * @param[in,out] mu Gravitational parameter.
    * @param[in] factors Structure containing the factors for dimensionalisation.
    */ 
    void cartesian_dimensionalise(double &t, Vector6 &RV, double &mu, const DimensionalFactors &factors);

}

#endif