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
     * @author Max Hallgarten La Casta
     * @date 2022-01-19
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
     * @author Max Hallgarten La Casta
     * @date 2022-01-19
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
     * @author Max Hallgarten La Casta
     * @date 2022-01-19
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
     * @author Max Hallgarten La Casta
     * @date 2022-01-19
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
     * @author Max Hallgarten La Casta
     * @date 2022-01-19
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
     * @author Max Hallgarten La Casta
     * @date 2022-01-19
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
     * @author Max Hallgarten La Casta
     * @date 2022-01-19
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
     * @author Max Hallgarten La Casta
     * @date 2022-01-19
     * 
     * @tparam T Numeric type.
     * @param[in] a First vector.
     * @param[in] b Second vector.
     * @return std::vector<T> Output vector.
     */
    template<class T>
    std::vector<T> cross3(const std::vector<T>& a, const std::vector<T>& b);

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Function to calculate the dot product of two vectors of polynomials with three elements.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-01-27
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] a First vector of polynomials.
     * @param[in] b Second vector of polynomials.
     * @return P<T> Dot product of the vectors of polynomials.
     */
    template<class T, template<class> class P>
    P<T> dot3(const std::vector<P<T>>& a, const std::vector<P<T>>& b);

    /**
     * @brief Function to calculate the length of a vector of polynomials with three elements.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-01-27
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] a Vector of polynomials.
     * @return P<T> Length of the vector of polynomials.
     */
    template<class T, template<class> class P>
    P<T> norm3(const std::vector<P<T>>& a);

    /**
     * @brief Calculate the cross product of two vectors of polynomials with three elements.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-01-31
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] a First vector of polynomials.
     * @param[in] b Second vector of polynomials.
     * @param[out] vecout Cross product of the vectors of polynomials.
     */
    template<class T, template<class> class P>
    void cross3(const std::vector<P<T>>& a, const std::vector<P<T>>& b, std::vector<P<T>>& vecout);

    /**
     * @brief Calculate the cross product of two vectors of polynomials with three elements.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-01-31
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] a First vector of polynomials.
     * @param[in] b Second vector of polynomials.
     * @return std::vector<P<T>> Cross product of the vectors of polynomials.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> cross3(const std::vector<P<T>>& a, const std::vector<P<T>>& b);

    #endif

}

#endif