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
#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../../include/vector/arithmeticoverloads.h"

namespace thames::vector::arithmeticoverloads {

    ////////////
    // Arrays //
    ////////////

    template<class T, std::size_t S>
    std::array<T, S> operator+(const std::array<T, S>& a, const std::array<T, S>& b){
        std::array<T, S> c;

        for(std::size_t ii=0; ii<S; ii++)
            c[ii] = a[ii] + b[ii];

        return c;
    }
    template std::array<double, 3> operator+<double, 3>(const std::array<double, 3>& a, const std::array<double, 3>& b);

    template<class T, std::size_t S>
    std::array<T, S> operator-(const std::array<T, S>& a, const std::array<T, S>& b){
        std::array<T, S> c;

        for(std::size_t ii=0; ii<S; ii++)
            c[ii] = a[ii] - b[ii];

        return c;
    }
    template std::array<double, 3> operator-<double, 3>(const std::array<double, 3>& a, const std::array<double, 3>& b);

    template<class T, std::size_t S>
    std::array<T, S> operator*(const T& a, const std::array<T, S>& b){
        std::array<T, S> c;

        for(std::size_t ii=0; ii<S; ii++)
            c[ii] = a*b[ii];

        return c;
    }
    template std::array<double, 3> operator*<double, 3>(const double& a, const std::array<double, 3>& b);

    template<class T, std::size_t S>
    std::array<T, S> operator*(const std::array<T, S>& b, const T& a){
        return a*b;
    }
    template std::array<double, 3> operator*<double, 3>(const std::array<double, 3>& b, const double& a);

    template<class T, std::size_t S>
    std::array<T, S> operator/(const std::array<T, S>& b, const T& a){
        std::array<T, S> c;

        for(std::size_t ii=0; ii<S; ii++)
            c[ii] = b[ii]/a;

        return c;
    }
    template std::array<double, 3> operator/<double, 3>(const std::array<double, 3>& b, const double& a);

    /////////////
    // Vectors //
    /////////////

    template<class T>
    std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b){
        std::vector<T> c(a.size());

        for(std::size_t ii=0; ii<a.size(); ii++)
            c[ii] = a[ii] + b[ii];

        return c;
    }
    template std::vector<double> operator+<double>(const std::vector<double>& a, const std::vector<double>& b);

    template<class T>
    std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b){
        std::vector<T> c(a.size());

        for(std::size_t ii=0; ii<a.size(); ii++)
            c[ii] = a[ii] - b[ii];

        return c;
    }
    template std::vector<double> operator-<double>(const std::vector<double>& a, const std::vector<double>& b);

    template<class T>
    std::vector<T> operator*(const T& a, const std::vector<T>& b){
        std::vector<T> c(b.size());

        for(std::size_t ii=0; ii<b.size(); ii++)
            c[ii] = a*b[ii];

        return c;
    }
    template std::vector<double> operator*<double>(const double& a, const std::vector<double>& b);

    template<class T>
    std::vector<T> operator*(const std::vector<T>& b, const T& a){
        return a*b;
    }
    template std::vector<double> operator*<double>(const std::vector<double>& b, const double& a);

    template<class T>
    std::vector<T> operator/(const std::vector<T>& b, const T& a){
        std::vector<T> c(b.size());

        for(std::size_t ii=0; ii<b.size(); ii++)
            c[ii] = b[ii]/a;

        return c;
    }
    template std::vector<double> operator/<double>(const std::vector<double>& b, const double& a);

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::polynomial;

    template<class T, template<class> class P>
    std::vector<P<T>> operator+(const std::vector<P<T>>& a, const std::vector<P<T>>& b){
        std::vector<P<T>> c(a);

        for(std::size_t ii=0; ii<a.size(); ii++)
            c[ii] = a[ii] + b[ii];

        return c;
    }
    template std::vector<taylor_polynomial<double>> operator+(const std::vector<taylor_polynomial<double>>& a, const std::vector<taylor_polynomial<double>>& b);
    template std::vector<chebyshev_polynomial<double>> operator+(const std::vector<chebyshev_polynomial<double>>& a, const std::vector<chebyshev_polynomial<double>>& b);

    template<class T, template<class> class P>
    std::vector<P<T>> operator-(const std::vector<P<T>>& a, const std::vector<P<T>>& b){
        std::vector<P<T>> c(a);

        for(std::size_t ii=0; ii<a.size(); ii++)
            c[ii] = a[ii] - b[ii];

        return c;
    }
    template std::vector<taylor_polynomial<double>> operator-(const std::vector<taylor_polynomial<double>>& a, const std::vector<taylor_polynomial<double>>& b);
    template std::vector<chebyshev_polynomial<double>> operator-(const std::vector<chebyshev_polynomial<double>>& a, const std::vector<chebyshev_polynomial<double>>& b);

    template<class T, template<class> class P>
    std::vector<P<T>> operator*(const T& a, const std::vector<P<T>>& b){
        std::vector<P<T>> c(b);

        for(std::size_t ii=0; ii<b.size(); ii++)
            c[ii] = a*b[ii];

        return c;
    }
    template std::vector<taylor_polynomial<double>> operator*(const double& a, const std::vector<taylor_polynomial<double>>& b);
    template std::vector<chebyshev_polynomial<double>> operator*(const double& a, const std::vector<chebyshev_polynomial<double>>& b);

    template<class T, template<class> class P>
    std::vector<P<T>> operator*(const std::vector<P<T>>& b, const T& a){
        return a*b;
    }
    template std::vector<taylor_polynomial<double>> operator*(const std::vector<taylor_polynomial<double>>& b, const double& a);
    template std::vector<chebyshev_polynomial<double>> operator*(const std::vector<chebyshev_polynomial<double>>& b, const double& a);

    template<class T, template<class> class P>
    std::vector<P<T>> operator*(const P<T>& a, const std::vector<P<T>>& b){
        std::vector<P<T>> c(b);

        for(std::size_t ii=0; ii<b.size(); ii++)
            c[ii] = a*b[ii];

        return c;
    }
    template std::vector<taylor_polynomial<double>> operator*(const taylor_polynomial<double>& a, const std::vector<taylor_polynomial<double>>& b);
    template std::vector<chebyshev_polynomial<double>> operator*(const chebyshev_polynomial<double>& a, const std::vector<chebyshev_polynomial<double>>& b);

    template<class T, template<class> class P>
    std::vector<P<T>> operator*(const std::vector<P<T>>& b, const P<T>& a){
        return a*b;
    }
    template std::vector<taylor_polynomial<double>> operator*(const std::vector<taylor_polynomial<double>>& b, const taylor_polynomial<double>& a);
    template std::vector<chebyshev_polynomial<double>> operator*(const std::vector<chebyshev_polynomial<double>>& b, const chebyshev_polynomial<double>& a);

    template<class T, template<class> class P>
    std::vector<P<T>> operator/(const std::vector<P<T>>& b, const T& a){
        std::vector<P<T>> c(b);

        for(std::size_t ii=0; ii<b.size(); ii++)
            c[ii] = b[ii]/a;

        return c;
    }
    template std::vector<taylor_polynomial<double>> operator/(const std::vector<taylor_polynomial<double>>& b, const double& a);
    template std::vector<chebyshev_polynomial<double>> operator/(const std::vector<chebyshev_polynomial<double>>& b, const double& a);

    template<class T, template<class> class P>
    std::vector<P<T>> operator/(const std::vector<P<T>>& b, const P<T>& a){
        std::vector<P<T>> c(b);

        for(std::size_t ii=0; ii<b.size(); ii++)
            c[ii] = b[ii]/a;

        return c;
    }
    template std::vector<taylor_polynomial<double>> operator/(const std::vector<taylor_polynomial<double>>& b, const taylor_polynomial<double>& a);
    template std::vector<chebyshev_polynomial<double>> operator/(const std::vector<chebyshev_polynomial<double>>& b, const chebyshev_polynomial<double>& a);

    #endif

}