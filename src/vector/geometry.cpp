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

#include <array>
#include <cmath>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../../include/vector/geometry.h"

namespace thames::vector::geometry{

    ///////////
    // Reals //
    ///////////

    template<class T>
    T dot3(const std::vector<T>& a, const std::vector<T>& b){
        // Return dot product
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    template double dot3<double>(const std::vector<double>&, const std::vector<double>&);

    template<class T>
    T norm3(const std::vector<T>& a){
        // Return square root of the dot product of the vector and itself
        return sqrt(dot3<T>(a, a));
    }
    template double norm3<double>(const std::vector<double>&);

    template<class T>
    void cross3(const std::vector<T>& a, const std::vector<T>& b, std::vector<T>& vecout){
        // Calculate cross product
        vecout[0] = a[1]*b[2] - a[2]*b[1];
        vecout[1] = a[2]*b[0] - a[0]*b[2];
        vecout[2] = a[0]*b[1] - a[1]*b[0];
    }
    template void cross3<double>(const std::vector<double>&, const std::vector<double>&, std::vector<double>&);

    template<class T>
    std::vector<T> cross3(const std::vector<T>& a, const std::vector<T>& b){
        // Declare output vector
        std::vector<T> vecout(3);

        // Calculate cross product
        cross3<T>(a, b, vecout);

        // Return cross product
        return vecout;
    }
    template std::vector<double> cross3<double>(const std::vector<double>&, const std::vector<double>&);

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::polynomial;

    template<class T, template<class> class P>
    P<T> dot3(const std::vector<P<T>>& a, const std::vector<P<T>>& b){
        // Return dot product
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    template taylor_polynomial<double> dot3(const std::vector<taylor_polynomial<double>>&, const std::vector<taylor_polynomial<double>>&);
    template chebyshev_polynomial<double> dot3(const std::vector<chebyshev_polynomial<double>>&, const std::vector<chebyshev_polynomial<double>>&);

    template<class T, template<class> class P>
    P<T> norm3(const std::vector<P<T>>& a){
        // Return square root of the dot product of the vector and itself
        return sqrt(dot3(a, a));
    }
    template taylor_polynomial<double> norm3(const std::vector<taylor_polynomial<double>>&);
    template chebyshev_polynomial<double> norm3(const std::vector<chebyshev_polynomial<double>>&);

    template<class T, template<class> class P>
    void cross3(const std::vector<P<T>>& a, const std::vector<P<T>>& b, std::vector<P<T>>& vecout){
        // Calculate cross product
        vecout[0] = a[1]*b[2] - a[2]*b[1];
        vecout[1] = a[2]*b[0] - a[0]*b[2];
        vecout[2] = a[0]*b[1] - a[1]*b[0];
    }
    template void cross3(const std::vector<taylor_polynomial<double>>&, const std::vector<taylor_polynomial<double>>&, std::vector<taylor_polynomial<double>>&);
    template void cross3(const std::vector<chebyshev_polynomial<double>>&, const std::vector<chebyshev_polynomial<double>>&, std::vector<chebyshev_polynomial<double>>&);

    template<class T, template<class> class P>
    std::vector<P<T>> cross3(const std::vector<P<T>>& a, const std::vector<P<T>>& b){
        // Declare output vector
        std::vector<P<T>> vecout(a);

        // Calculate cross product
        cross3(a, b, vecout);

        // Return cross product
        return vecout;
    }
    template std::vector<taylor_polynomial<double>> cross3(const std::vector<taylor_polynomial<double>>&, const std::vector<taylor_polynomial<double>>&);
    template std::vector<chebyshev_polynomial<double>> cross3(const std::vector<chebyshev_polynomial<double>>&, const std::vector<chebyshev_polynomial<double>>&);

    #endif

}