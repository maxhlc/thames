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
     * @author Max Hallgarten La Casta
     * @date 2022-04-29
     * 
     * @tparam T Numeric type.
     * @param[in] RV Dimensional Cartesian state vector.
     * @param[in] factors Structure containing the factors for non-dimensionalisation.
     * @return std::array<T, 6> Non-dimensional Cartesian state vector.
     */
    template<class T>
    std::array<T, 6> cartesian_nondimensionalise(const std::array<T, 6>& RV, const DimensionalFactors<T>& factors);

    /**
     * @brief Dimensionalise Cartesian state.
     *
     * @author Max Hallgarten La Casta
     * @date 2022-04-29
     * 
     * @tparam T Numeric type.
     * @param[in] RVnd Non-dimensional Cartesian state vector.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     * @return std::array<T, 6> Dimensional Cartesian state vector.
     */ 
    template<class T>
    std::array<T, 6> cartesian_dimensionalise(const std::array<T, 6>& RVnd, const DimensionalFactors<T>& factors);

    /**
     * @brief Non-dimensionalise GEqOE state.
     *
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @param[in] geqoe Dimensional GEqOE state vector.
     * @param[in] factors Structure containing the factors for non-dimensionalisation.
     * @return std::array<T, 6> Non-dimensional GEqOE state vector.
     */
    template<class T>
    std::array<T, 6> geqoe_nondimensionalise(const std::array<T, 6>& geqoe, const DimensionalFactors<T>& factors);

    /**
     * @brief Dimensionalise GEqOE state.
     *
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @param[in] geqoend Non-dimensional GEqOE state vector.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     * @return std::array<T, 6> Dimensional GEqOE state vector.
     */ 
    template<class T>
    std::array<T, 6> geqoe_dimensionalise(const std::array<T, 6>& geqoend, const DimensionalFactors<T>& factors);

    /**
     * @brief Calculate non-dimensionalisation factors.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-03
     * 
     * @tparam T Numeric type.
     * @param[in] RV Cartesian state vector. 
     * @param[in] mu Gravitational parameter.
     * @return DimensionalFactors<T> Non-dimensionalisation factors.
     */
    template<class T>
    DimensionalFactors<T> calculate_factors(const std::array<T, 6>& RV, const T& mu);

    /////////////
    // Vectors //
    /////////////

    /**
     * @brief Non-dimensionalise Cartesian state.
     *
     * @author Max Hallgarten La Casta
     * @date 2022-04-29
     * 
     * @tparam T Numeric type.
     * @param[in] RV Dimensional Cartesian state vector.
     * @param[in] factors Structure containing the factors for non-dimensionalisation.
     * @return std::vector<T> Non-dimensionalised Cartesian state vector.
     */
    template<class T>
    std::vector<T> cartesian_nondimensionalise(const std::vector<T>& RV, const DimensionalFactors<T>& factors);

    /**
     * @brief Dimensionalise Cartesian state.
     *
     * @author Max Hallgarten La Casta
     * @date 2022-04-29
     * 
     * @tparam T Numeric type.
     * @param[in] RVnd Non-dimensional Cartesian state vector.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     * @return std::vector<T> Dimensional Cartesian state vector.
     */ 
    template<class T>
    std::vector<T> cartesian_dimensionalise(const std::vector<T>& RVnd, const DimensionalFactors<T>& factors);

    /**
     * @brief Non-dimensionalise GEqOE state.
     *
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @param[in] geqoe Dimensional GEqOE state vector.
     * @param[in] factors Structure containing the factors for non-dimensionalisation.
     * @return std::array<T, 6> Non-dimensional GEqOE state vector.
     */
    template<class T>
    std::vector<T> geqoe_nondimensionalise(const std::vector<T>& geqoe, const DimensionalFactors<T>& factors);

    /**
     * @brief Dimensionalise GEqOE state.
     *
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @param[in] geqoend Non-dimensional GEqOE state vector.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     * @return std::array<T, 6> Dimensional GEqOE state vector.
     */ 
    template<class T>
    std::vector<T> geqoe_dimensionalise(const std::vector<T>& geqoend, const DimensionalFactors<T>& factors);

    /**
     * @brief Calculate non-dimensionalisation factors.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-04-29
     * 
     * @tparam T Numeric type.
     * @param[in] RV Cartesian state vector. 
     * @param[in] mu Gravitational parameter.
     * @return DimensionalFactors<T> Non-dimensionalisation factors.
     */
    template<class T>
    DimensionalFactors<T> calculate_factors(const std::vector<T>& RV, const T& mu);

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Non-dimensionalise Cartesian state polynomial.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-04-29
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] RV Dimensional Cartesian state vector.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     * @return std::vector<P<T>> Non-dimensional Cartesian state vector.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> cartesian_nondimensionalise(const std::vector<P<T>>& RV, const DimensionalFactors<T>& factors);

    /**
     * @brief Dimensionalise Cartesian state polynomial.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-04-29
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] RVnd Non-dimensional Cartesian state vector.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     * @return std::vector<P<T>> Dimensional Cartesian state vector.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> cartesian_dimensionalise(const std::vector<P<T>>& RVnd, const DimensionalFactors<T>& factors);

    /**
     * @brief Non-dimensionalise GEqOE state polynomial.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] geqoe Dimensional GEqOE state vector.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     * @return std::vector<P<T>> Non-dimensional GEqOE state vector.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> geqoe_nondimensionalise(const std::vector<P<T>>& geqoe, const DimensionalFactors<T>& factors);

    /**
     * @brief Dimensionalise GEqOE state polynomial.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] RVnd Non-dimensional GEqOE state vector.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     * @return std::vector<P<T>> Dimensional GEqOE state vector.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> geqoe_dimensionalise(const std::vector<P<T>>& geqoend, const DimensionalFactors<T>& factors);   

    /**
     * @brief Calculate non-dimensionalisation factors.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-04-29
     * 
     * @tparam T Numeric type.
     * @param[in] RV Cartesian state vector. 
     * @param[in] mu Gravitational parameter.
     * @return DimensionalFactors<T> Non-dimensionalisation factors.
     */
    template<class T, template<class> class P>
    DimensionalFactors<T> calculate_factors(const std::vector<P<T>>& RV, const T& mu);

    #endif

}

#endif