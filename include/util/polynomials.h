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

#ifndef THAMES_UTIL_POLYNOMIALS
#define THAMES_UTIL_POLYNOMIALS

#include <vector>

namespace thames::util::polynomials {

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Evaluate vector of polynomials at a point.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-02-09
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] polynomials Vector of polynomials.
     * @param[in] x Evaluation point. 
     * @return std::vector<T> Point state vector.
     */
    template<class T, template<class> class P>
    std::vector<T> evaluate_polynomials(const std::vector<P<T>>& polynomials, const std::vector<T>& x);

    /**
     * @brief Evaluate vector of polynomials at a set of points.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-02-18
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] polynomials Vector of polynomials.
     * @param[in] x Evaluation points. 
     * @return std::vector<std::vector<T>> Point state vectors.
     */   
    template<class T, template<class> class P>
    std::vector<std::vector<T>> evaluate_polynomials(const std::vector<P<T>>& polynomials, const std::vector<std::vector<T>>& x);

    #endif

}

#endif