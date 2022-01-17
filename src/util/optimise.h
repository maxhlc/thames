#ifndef THAMES_UTIL_OPTIMISE
#define THAMES_UTIL_OPTIMISE

#include <functional>

namespace thames::util::optimise{

    /**
     * @brief Golden section search to find local minima.
     * 
     * @details Method as described by: https://en.wikipedia.org/wiki/Golden-section_search#Iterative_algorithm
     * 
     * @tparam T Numeric type.
     * @param[in] func Scalar function to be minimised.
     * @param[in] a Left hand boundary.
     * @param[in] b Right hand boundary.
     * @param[in] tol Solver tolerance.
     * @return T Argument of the minimum.
     */
    template<class T>
    T golden_section_search(std::function<T (T)> func, T a, T b, T tol = 1e-10);

}

#endif