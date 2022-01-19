#ifndef THAMES_VECTOR_GEOMETRY
#define THAMES_VECTOR_GEOMETRY

#include <array>
#include <vector>

namespace thames::vector::geometry{

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
    T norm3(const std::array<T, 3>& a);

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
     * @brief Function to return the slice of a vector.
     * 
     * @tparam T Numeric type.
     * @tparam Si Input size.
     * @tparam So Output size.
     * @param[in] v Input vector for slicing.
     * @param[in] a Lower index.
     * @param[in] b Upper index.
     * @return std::array<T, So> Output vector which has been sliced.
     */
    template<class T, const std::size_t Si, const std::size_t So>
    std::array<T, So> slice(const std::array<T, Si>& v, const unsigned int a, const unsigned int b);

    /**
     * @brief Function to return the slice of a vector.
     * 
     * @tparam T Numeric type.
     * @param[in] v Input vector for slicing.
     * @param[in] a Lower index.
     * @param[in] b Upper index.
     * @return std::vector<T> Output vector which has been sliced.
     */
    template<class T>
    std::vector<T> slice(const std::vector<T>& v, const unsigned int a, const unsigned int b);

    template<class T>
    std::array<T, 3> mult3(const T a, const std::array<T, 3>& vec);

    template<class T>
    std::vector<T> mult3(const T a, const std::vector<T>& vec);

    template<class T>
    std::array<T, 3> cross3(const std::array<T, 3>& a, const std::array<T, 3>& b);

    template<class T>
    std::vector<T> cross3(const std::vector<T>& a, const std::vector<T>& b);

}

#endif