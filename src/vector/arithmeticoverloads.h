#ifndef THAMES_VECTOR_ARITHMETICOVERLOADS
#define THAMES_VECTOR_ARITHMETICOVERLOADS

#include <array>
#include <vector>

namespace thames::vector::arithmeticoverloads {

    ////////////
    // Arrays //
    ////////////

    /**
     * @brief Element-wise addition of two arrays.
     * 
     * @tparam T Numeric type.
     * @tparam S Array size.
     * @param[in] a First array.
     * @param[in] b Second array.
     * @return std::array<T, S> Output array.
     */
    template<class T, std::size_t S>
    std::array<T, S> operator+(const std::array<T, S>& a, const std::array<T, S>& b);

    /**
     * @brief Element-wise subtraction of two arrays.
     * 
     * @tparam T Numeric type.
     * @tparam S Array size.
     * @param[in] a First array.
     * @param[in] b Second array.
     * @return std::array<T, S> Output array.
     */
    template<class T, std::size_t S>
    std::array<T, S> operator-(const std::array<T, S>& a, const std::array<T, S>& b);

    /**
     * @brief Scalar multiplication of an array.
     * 
     * @tparam T Numeric type.
     * @tparam S Array size.
     * @param[in] a Scalar.
     * @param[in] b Array.
     * @return std::array<T, S> Output array.
     */
    template<class T, std::size_t S>
    std::array<T, S> operator*(const T& a, const std::array<T, S>& b);

    /**
     * @brief Scalar multiplication of an array.
     * 
     * @tparam T Numeric type.
     * @tparam S Array size.
     * @param[in] b Array.
     * @param[in] a Scalar.
     * @return std::array<T, S> Output array.
     */
    template<class T, std::size_t S>
    std::array<T, S> operator*(const std::array<T, S>& b, const T& a);

    /**
     * @brief Scalar division of an array.
     * 
     * @tparam T Numeric type.
     * @tparam S Array size.
     * @param[in] b Array.
     * @param[in] a Scalar.
     * @return std::array<T, S> Output array.
     */
    template<class T, std::size_t S>
    std::array<T, S> operator/(const std::array<T, S>& b, const T& a);

    /////////////
    // Vectors //
    /////////////

    /**
     * @brief Element-wise addition of two vectors.
     * 
     * @tparam T Numeric type.
     * @param[in] a First vector.
     * @param[in] b Second vector.
     * @return std::vector<T> Output vector.
     */
    template<class T>
    std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b);

    /**
     * @brief Element-wise subtraction of two vectors.
     * 
     * @tparam T Numeric type.
     * @param[in] a First vector.
     * @param[in] b Second vector.
     * @return std::vector<T> Output vector.
     */
    template<class T>
    std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b);

    /**
     * @brief Scalar multiplication of a vector.
     * 
     * @tparam T Numeric type.
     * @param[in] a Scalar.
     * @param[in] b Array.
     * @return std::vector<T> Output array.
     */
    template<class T>
    std::vector<T> operator*(const T& a, const std::vector<T>& b);

    /**
     * @brief Scalar multiplication of a vector.
     * 
     * @tparam T Numeric type.
     * @param[in] b Array.
     * @param[in] a Scalar.
     * @return std::vector<T> Output array.
     */
    template<class T>
    std::vector<T> operator*(const std::vector<T>& b, const T& a);

    /**
     * @brief Scalar division of a vector.
     * 
     * @tparam T Numeric type.
     * @param[in] b Array.
     * @param[in] a Scalar.
     * @return std::vector<T> Output vector.
     */
    template<class T>
    std::vector<T> operator/(const std::vector<T>& b, const T& a);

}

#endif