#ifndef THAMES_CONVERSIONS_DIMENSIONAL
#define THAMES_CONVERSIONS_DIMENSIONAL

#include "state.h"
#include "../types.h"

using namespace thames::types;

namespace thames::conversions::dimensional{

    /**
    * @brief Non-dimensionalise Cartesian state.
    *
    * @param t Current physical time.
    * @param RV Current Cartesian state vector (position and velocity).
    * @param mu Gravitational parameter.
    * @param factors Structure containing the factors for non-dimensionalisation.
    */ 
    void cartesian_nondimensionalise(double &t, Vector6 &RV, double &mu, DimensionalFactors &factors);

    /**
    * @brief Dimensionalise Cartesian state.
    *
    * @param t Current physical time.
    * @param RV Current Cartesian state vector (position and velocity).
    * @param mu Gravitational parameter.
    * @param factors Structure containing the factors for dimensionalisation.
    */ 
    void cartesian_dimensionalise(double &t, Vector6 &RV, double &mu, const DimensionalFactors &factors);

}

#endif