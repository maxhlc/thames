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

#include <cmath>
#include <functional>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../../include/util/optimise.h"
#include "../../include/util/root.h"

namespace thames::util::root{

    ///////////
    // Reals //
    ///////////

    template<class T>
    T golden_section_search(std::function<T (T)> func, T a, T b, T tol){
        // Declare function for minimisation (root at minimum of absolute of the function)
        std::function<T (T)> f = [func](T x) {return fabs(func(x));};

        // Minimise function using the golden section search
        T x0 = thames::util::optimise::golden_section_search(f, a, b, tol);

        // Return root
        return x0;
    }
    template double golden_section_search<double>(std::function<double (double)>, double, double, double);

    template<class T>
    T newton_raphson(const std::function<T (T)>& func, const std::function<T (T)>& dfunc, T xn, T tol){
        // Declare approximation variable
        T xn1;

        // Set converged flag to false
        bool converged = false;

        // Iterate until converged
        while(!converged){
            // Update approximation
            xn1 = xn - func(xn)/dfunc(xn);

            // Converged if update is smaller than tolerance
            if(fabs(xn1 - xn) < tol)
                converged = true;

            // Update previous approximation
            xn = xn1;
        }

        // Return root
        return xn1;
    }
    template double newton_raphson<double>(const std::function<double (double)>&, const std::function<double (double)>&, double, double);

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::polynomial;

    template<class T, template<class> class P>
    P<T> newton_raphson(const std::function<P<T> (P<T>)>& func, const std::function<P<T> (P<T>)>& dfunc, P<T> xn, T tol){
        // Declare approximation variable
        P<T> xn1(xn);

        // Set converged flag to false
        bool converged = false;

        // Iterate until converged
        while(!converged){
            // Update approximation
            xn1 = xn - func(xn)/dfunc(xn);

            // Converged if update is smaller than tolerance
            if(fabs((xn1 - xn).get_coeffs()[0]) < tol)
                converged = true;

            // Update previous approximation
            xn = xn1;
        }

        // Return root
        return xn1;
    }
    template taylor_polynomial<double> newton_raphson(const std::function<taylor_polynomial<double> (taylor_polynomial<double>)>&, const std::function<taylor_polynomial<double> (taylor_polynomial<double>)>&, taylor_polynomial<double>, double);
    template chebyshev_polynomial<double> newton_raphson(const std::function<chebyshev_polynomial<double> (chebyshev_polynomial<double>)>&, const std::function<chebyshev_polynomial<double> (chebyshev_polynomial<double>)>&, chebyshev_polynomial<double>, double);

    #endif

}