/*
MIT License

Copyright (c) 2021-2022 Max Hallgarten La Casta

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

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