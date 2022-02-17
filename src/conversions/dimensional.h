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

#ifndef THAMES_CONVERSIONS_DIMENSIONAL
#define THAMES_CONVERSIONS_DIMENSIONAL

#include <array>
#include <vector>

namespace thames::conversions::dimensional{

    /**
     * @brief Structure to contain factors for (non)dimensionalisation of Cartesian states.
     * 
     * @tparam T Numeric type.
     */
    template <class T>
    struct DimensionalFactors {
        T time = 1.0;
        T length = 1.0;
        T velocity = 1.0;
        T grav = 1.0;
    };

    ////////////
    // Arrays //
    ////////////

    /**
     * @brief Non-dimensionalise Cartesian state.
     *
     * @tparam T Numeric type.
     * @param[in,out] t Current physical time.
     * @param[in,out] RV Current Cartesian state vector (position and velocity).
     * @param[in,out] mu Gravitational parameter.
     * @param[out] factors Structure containing the factors for non-dimensionalisation.
     */
    template<class T>
    void cartesian_nondimensionalise(T& t, std::array<T, 6>& RV, T& mu, DimensionalFactors<T>& factors);

    /**
     * @brief Dimensionalise Cartesian state.
     *
     * @tparam T Numeric type.
     * @param[in,out] t Current physical time.
     * @param[in,out] RV Current Cartesian state vector (position and velocity).
     * @param[in,out] mu Gravitational parameter.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     */ 
    template<class T>
    void cartesian_dimensionalise(T& t, std::array<T, 6>& RV, T& mu, const DimensionalFactors<T>& factors);

    /**
     * @brief Calculate non-dimensionalisation factors.
     * 
     * @tparam T Numeric type.
     * @param[in] RV Cartesian state vector. 
     * @param[in] mu Gravitational parameter.
     * @return DimensionalFactors<T> Non-dimensionalisation factors.
     */
    template<class T>
    DimensionalFactors<T> calculate_factors(const std::array<T, 6> RV, const T mu);

    /////////////
    // Vectors //
    /////////////

    /**
     * @brief Non-dimensionalise Cartesian state.
     *
     * @tparam T Numeric type.
     * @param[in,out] t Current physical time.
     * @param[in,out] RV Current Cartesian state vector (position and velocity).
     * @param[in,out] mu Gravitational parameter.
     * @param[out] factors Structure containing the factors for non-dimensionalisation.
     */
    template<class T>
    void cartesian_nondimensionalise(T& t, std::vector<T>& RV, T& mu, DimensionalFactors<T>& factors);

    /**
     * @brief Dimensionalise Cartesian state.
     *
     * @tparam T Numeric type.
     * @param[in,out] t Current physical time.
     * @param[in,out] RV Current Cartesian state vector (position and velocity).
     * @param[in,out] mu Gravitational parameter.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     */ 
    template<class T>
    void cartesian_dimensionalise(T& t, std::vector<T>& RV, T& mu, const DimensionalFactors<T>& factors);

    /**
     * @brief Calculate non-dimensionalisation factors.
     * 
     * @tparam T Numeric type.
     * @param[in] RV Cartesian state vector. 
     * @param[in] mu Gravitational parameter.
     * @return DimensionalFactors<T> Non-dimensionalisation factors.
     */
    template<class T>
    DimensionalFactors<T> calculate_factors(const std::vector<T> RV, const T mu);

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Non-dimensionalise Cartesian state polynomial.
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in,out] t Current physical time.
     * @param[in,out] RV Current Cartesian state vector (position and velocity).
     * @param[in,out] mu Gravitational parameter.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     */
    template<class T, template<class> class P>
    void cartesian_nondimensionalise(T& t, std::vector<P<T>>& RV, T& mu, DimensionalFactors<T>& factors);

    /**
     * @brief Dimensionalise Cartesian state polynomial.
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in,out] t Current physical time.
     * @param[in,out] RV Current Cartesian state vector (position and velocity).
     * @param[in,out] mu Gravitational parameter.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     */
    template<class T, template<class> class P>
    void cartesian_dimensionalise(T& t, std::vector<P<T>>& RV, T& mu, DimensionalFactors<T>& factors);

    #endif

}

#endif