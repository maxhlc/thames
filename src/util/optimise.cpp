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

#include "../../include/util/optimise.h"

namespace thames::util::optimise{

    template<class T>
    T golden_section_search(std::function<T (T)> func, T a, T b, T tol){
        // Set parameters
        const T gr = 0.5*(sqrt(5.0) + 1);

        // Calculate test points
        T c = b - (b - a)/gr;
        T d = a + (b - a)/gr;

        // Iterate until convergence
        while (fabs(b - a) > tol) {
            // Choose test points based on function values
            if (func(c) < func(d)) {
                b = d;
            } else {
                a = c;
            }

            // Recompute test points
            c = b - (b - a)/gr;
            d = a + (b - a)/gr;           
        }

        // Return converged value
        return (0.5*(a + b));
    }
    template double golden_section_search<double>(std::function<double (double)>, double, double, double);

}