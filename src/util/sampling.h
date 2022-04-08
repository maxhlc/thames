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

#ifndef THAMES_UTIL_SAMPLING
#define THAMES_UTIL_SAMPLING

namespace thames::util::sampling {

    /**
     * @brief Calculate evenly spaced points over the [a,b] interval.
     * 
     * Returns the mid-point if one point is requested.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-02-01
     * 
     * @tparam T Numeric type.
     * @param[in] a Start value.
     * @param[in] b End value.
     * @param[in] n Number of points.
     * @return std::vector<T> Evenly space points.
     */
    template<class T>
    std::vector<T> linspace(const T a, const T b, const unsigned int n);

    /**
     * @brief Generate permuations of Cartesian states from sets of points in each state variable.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-02-01
     * 
     * @tparam T Numeric type.
     * @param[in] points Vector of vectors of ranges for each Cartesian state (i.e. X, Y, Z, VX, VY, VZ)
     * @return std::vector<std::vector<T>> 
     */
    template<class T>
    std::vector<std::vector<T>> cartesian_permutations(const std::vector<std::vector<T>>& points);

}

#endif