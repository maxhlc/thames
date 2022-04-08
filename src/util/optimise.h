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

#ifndef THAMES_UTIL_OPTIMISE
#define THAMES_UTIL_OPTIMISE

#include <functional>

namespace thames::util::optimise{

    /**
     * @brief Golden section search to find local minima.
     * 
     * @details Method as described by: https://en.wikipedia.org/wiki/Golden-section_search#Iterative_algorithm
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-01-17
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