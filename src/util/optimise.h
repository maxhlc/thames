#ifndef THAMES_UTIL_OPTIMISE
#define THAMES_UTIL_OPTIMISE

#include <functional>
#include <math.h>

namespace thames::util::optimise{

    /**
     * @brief Golden section search to find local minima.
     * 
     * @details Method as described by: https://en.wikipedia.org/wiki/Golden-section_search#Iterative_algorithm
     * 
     * @tparam real Type for a real number.
     * @param[in] func Scalar function to be minimised.
     * @param[in] a Left hand boundary.
     * @param[in] b Right hand boundary.
     * @param[in] tol Solver tolerance.
     * @return real Argument of the minimum.
     */
    template<class real>
    real golden_section_search(std::function<real (real)> func, real a, real b, real tol = 1e-10);

}

#endif