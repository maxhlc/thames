#ifndef THAMES_PROPAGATORS_GEQOE
#define THAMES_PROPAGATORS_GEQOE

#include "../types.h"

using namespace thames::types;

namespace thames::propagators::geqoe{

    void derivative(const Vector6 &geqoe, Vector6 &geqoedot, const double t, const double &mu, const Potential &U_func, const PotentialDerivative &Ut_func, const Force &F_func, const Force &P_func);

}

#endif