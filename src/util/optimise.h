#ifndef THAMES_UTIL_OPTIMISE
#define THAMES_UTIL_OPTIMISE

#include <functional>
#include <math.h>

namespace thames::util::optimise{

    /**
     * @brief Golden section search to find local minima.
     * 
     * @details Method descrived by: https://en.wikipedia.org/wiki/Golden-section_search#Iterative_algorithm
     * 
     * @param[in] func Scalar function to be minimised.
     * @param[in] a Left hand boundary.
     * @param[in] b Right hand boundary.
     * @return double Argument of the minimum.
     */
    double golden_section_search(std::function<double (double)> func, double a, double b);

}

#endif