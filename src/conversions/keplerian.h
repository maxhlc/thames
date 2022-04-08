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

#ifndef THAMES_CONVERSIONS_KEPLERIAN
#define THAMES_CONVERSIONS_KEPLERIAN

#include <array>
#include <vector>

namespace thames::conversions::keplerian{

    ////////////
    // Arrays //
    ////////////

    /**
     * @brief Convert from Cartesian state to traditional Keplerian elements.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-01-24
     * 
     * @tparam T Numeric type.
     * @param[in] RV Cartesian state.
     * @param[in] mu Gravitational parameter.
     * @return std::array<T, 6> Keplerian elements state.
     */
    template<class T>
    std::array<T, 6> cartesian_to_keplerian(const std::array<T, 6>& RV, const T& mu);

    /**
     * @brief Convert from traditional Keplerian elements to Cartesian state.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-01-24
     * 
     * @tparam T Numeric type.
     * @param[in] keplerian Keplerian elements state.
     * @param[in] mu Gravitational parameter.
     * @return std::array<T, 6> Cartesian state.
     */
    template<class T>
    std::array<T, 6> keplerian_to_cartesian(const std::array<T, 6>& keplerian, const T& mu);

    /////////////
    // Vectors //
    /////////////

    /**
     * @brief Convert from Cartesian state to traditional Keplerian elements.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-01-24
     * 
     * @tparam T Numeric type.
     * @param[in] RV Cartesian state.
     * @param[in] mu Gravitational parameter.
     * @return std::vector<T> Keplerian elements state.
     */
    template<class T>
    std::vector<T> cartesian_to_keplerian(const std::vector<T>& RV, const T& mu);

    /**
     * @brief Convert from traditional Keplerian elements to Cartesian state.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-01-24
     * 
     * @tparam T Numeric type.
     * @param[in] keplerian Keplerian elements state.
     * @param[in] mu Gravitational parameter.
     * @return std::vector<T> Cartesian state.
     */
    template<class T>
    std::vector<T> keplerian_to_cartesian(const std::vector<T>& keplerian, const T& mu);

}

#endif