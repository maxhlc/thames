#ifndef THAMES_UTIL_VECTOR
#define THAMES_UTIL_VECTOR

#include "../types.h"

using namespace thames::types;

namespace thames::util::vector{

    template<class real, class vector>
    real dot3(vector a, vector b);

    template<class real, class vector>
    real norm3(vector a);

    template<class vectorout, class vectorin, class integer>
    vectorout slice(vectorin v, integer a, integer b);

}

#endif