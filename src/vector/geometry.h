#ifndef THAMES_VECTOR_GEOMETRY
#define THAMES_VECTOR_GEOMETRY

#include <array>
#include <vector>

namespace thames::vector::geometry{

    ////////////
    // Arrays //
    ////////////

    /**
     * @brief Function to calculate the dot product of two vectors with three elements.
     * 
     * @tparam T Numeric type.
     * @param[in] a First vector.
     * @param[in] b Second vector.
     * @return T Dot product of the vectors.
     */
    template<class T>
    T dot3(const std::array<T, 3>& a, const std::array<T, 3>& b);

    /**
     * @brief Function to calculate the length of a vector with three elements.
     * 
     * @tparam T Numeric type.
     * @param[in] a Vector.
     * @return T Length of the vector.
     */
    template<class T>
    T norm3(const std::array<T, 3>& a);

    /**
     * @brief Calculate the cross product of two vectors with three elements.
     * 
     * @tparam T Numeric type.
     * @param[in] a First vector.
     * @param[in] b Second vector.
     * @param[out] vecout Output vector.
     */
    template<class T>
    void cross3(const std::array<T, 3>& a, const std::array<T, 3>& b, std::array<T,3>& vecout);

    /**
     * @brief Calculate the cross product of two vectors with three elements.
     * 
     * @tparam T Numeric type.
     * @param[in] a First vector.
     * @param[in] b Second vector.
     * @return std::array<T, 3> Output vector.
     */
    template<class T>
    std::array<T, 3> cross3(const std::array<T, 3>& a, const std::array<T, 3>& b);

    /////////////
    // Vectors //
    /////////////

    /**
     * @brief Function to calculate the dot product of two vectors with three elements.
     * 
     * @tparam T Numeric type.
     * @param[in] a First vector.
     * @param[in] b Second vector.
     * @return T Dot product of the vectors.
     */
    template<class T>
    T dot3(const std::vector<T>& a, const std::vector<T>& b);

    /**
     * @brief Function to calculate the length of a vector with three elements.
     * 
     * @tparam T Numeric type.
     * @param[in] a Vector.
     * @return T Length of the vector.
     */
    template<class T>
    T norm3(const std::vector<T>& a);

    /**
     * @brief Calculate the cross product of two vectors with three elements.
     * 
     * @tparam T Numeric type.
     * @param[in] a First vector.
     * @param[in] b Second vector.
     * @param[out] vecout Output vector.
     */
    template<class T>
    void cross3(const std::vector<T>& a, const std::vector<T>& b, std::vector<T>& vecout);

    /**
     * @brief Calculate the cross product of two vectors with three elements.
     * 
     * @tparam T Numeric type.
     * @param[in] a First vector.
     * @param[in] b Second vector.
     * @return std::vector<T> Output vector.
     */
    template<class T>
    std::vector<T> cross3(const std::vector<T>& a, const std::vector<T>& b);

}

#endif