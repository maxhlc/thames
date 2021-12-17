#ifndef THAMES_UTIL_ROOT
#define THAMES_UTIL_ROOT

#include <functional>
#include <math.h>

namespace thames::util::root{

    /**
     * @brief Golden section search for root finding.
     * 
     * @param func Scalar function for root finding.
     * @param a Left hand boundary.
     * @param b Right hand boundary.
     * @return double Argument of the root.
     */
    double golden_section_search(std::function<double (double)> func, double a, double b);

}

#endif