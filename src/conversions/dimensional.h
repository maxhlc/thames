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

}

#endif