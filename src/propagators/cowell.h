#ifndef THAMES_PROPAGATORS_COWELL
#define THAMES_PROPAGATORS_COWELL

#include "../types.h"

using namespace thames::types;

namespace thames::propagators::cowell{

    void derivative(const Vector6 &RV, Vector6 &RVdot, const double t, const double &mu, const Force &F_func);

    Vector6 propagate(double tstart, double tend, double tstep, Vector6 RV, double mu, Force F_func);

}

#endif