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
     * @tparam T Numeric type.
     * @param[in] keplerian Keplerian elements state.
     * @param[in] mu Gravitational parameter.
     * @return std::vector<T> Cartesian state.
     */
    template<class T>
    std::vector<T> keplerian_to_cartesian(const std::vector<T>& keplerian, const T& mu);

}

#endif