#ifndef THAMES_UTIL_ROOT
#define THAMES_UTIL_ROOT

#include <functional>

namespace thames::util::root{

    ///////////
    // Reals //
    ///////////

    /**
     * @brief Golden section search for root finding.
     * 
     * @tparam T Numeric type.
     * @param[in] func Scalar function for root finding.
     * @param[in] a Left hand boundary.
     * @param[in] b Right hand boundary.
     * @param[in] tol Solver tolerance.
     * @return T Argument of the root.
     */
    template<class T>
    T golden_section_search(std::function<T (T)> func, T a, T b, T tol = 1e-10);

    /**
     * @brief Newton-Raphson method for root finding.
     * 
     * @tparam T Numeric type.
     * @param[in] func Scalar function for root finding. 
     * @param[in] dfunc Derivative of scalar function for root finding.
     * @param[in] xn Initial guess.
     * @param[in] tol Solver tolerance.
     * @return T Argument of the root.
     */
    template<class T>
    T newton_raphson(const std::function<T (T)> &func, const std::function<T (T)> &dfunc, T xn, T tol = 1e-10);

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Newton-Raphson method for root finding.
     * 
     * Convergence is calculated using the constant (index of zero) term.
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] func Scalar function for root finding. 
     * @param[in] dfunc Derivative of scalar function for root finding.
     * @param[in] xn Initial guess.
     * @param[in] tol Solver tolerance.
     * @return P<T> Argument of the root.
     */
    template<class T, template<class> class P>
    P<T> newton_raphson(const std::function<P<T> (P<T>)>& func, const std::function<P<T> (P<T>)>& dfunc, P<T> xn, T tol = 1e-10);

    #endif
}

#endif