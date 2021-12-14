#ifndef CONVERSIONS_UTIL
#define CONVERSIONS_UTIL

#include "../types.h"

using namespace thames::types;

namespace thames::conversions::util {

    Matrix33 rot_x(const double &ang);

    Matrix33 rot_y(const double &ang);
    
    Matrix33 rot_z(const double &ang);

}

#endif