#ifndef CONVERSIONS_DIMENSIONAL
#define CONVERSIONS_DIMENSIONAL

#include "state.h"
#include "../types.h"

using namespace thames::types;

namespace thames::conversions::dimensional{

    void cartesian_nondimensionalise(double &t, Vector6 &RV, double &mu, DimensionalFactors &factors);

    void cartesian_dimensionalise(double &t, Vector6 &RV, double &mu, const DimensionalFactors &factors);

}

#endif