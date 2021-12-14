#ifndef THAMES_PROPAGATORS_GEQOE
#define THAMES_PROPAGATORS_GEQOE

#include "../types.h"

using namespace thames::types;

namespace thames::propagators::geqoe{

    Vector6 derivative(const Vector6 &Y, Vector6 &Ydot, const double t);

}

#endif