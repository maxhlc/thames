#ifndef THAMES_UTIL_ROOT
#define THAMES_UTIL_ROOT

#include <functional>
#include <math.h>

namespace thames::util::root{

    /**
     * @brief Golden section search for root finding.
     * 
     * @tparam real Type for a real number.
     * @param[in] func Scalar function for root finding.
     * @param[in] a Left hand boundary.
     * @param[in] b Right hand boundary.
     * @param[in] tol Solver tolerance.
     * @return real Argument of the root.
     */
    template<class real>
    real golden_section_search(std::function<real (real)> func, real a, real b, real tol = 1e-10);

    /**
     * @brief Netwon-Raphson method for root finding.
     * 
     * @tparam real Type for a real number.
     * @param[in] func Scalar function for root finding. 
     * @param[in] dfunc Derivative of scalar function for root finding.
     * @param[in] xn Initial guess.
     * @param[in] tol Solver tolerance.
     * @return real Argument of the root.
     */
    template<class real>
    real newton_raphson(const std::function<real (real)> &func, const std::function<real (real)> &dfunc, real xn, real tol = 1e-10);

}

#endif