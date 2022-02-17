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

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Element-wise addition of two vectors of polynomials.
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] a First vector.
     * @param[in] b Second vector.
     * @return std::vector<P<T>> Output vector.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> operator+(const std::vector<P<T>>& a, const std::vector<P<T>>& b);

    /**
     * @brief Element-wise subtraction of two vectors of polynomials.
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] a First vector of polynomials.
     * @param[in] b Second vector of polynomials.
     * @return std::vector<P<T>> Output vector of polynomials.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> operator-(const std::vector<P<T>>& a, const std::vector<P<T>>& b);

    /**
     * @brief Multiplication of a vector of polynomials with a scalar.
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] a Scalar.
     * @param[in] b Vector of polynomials.
     * @return std::vector<P<T>> Output vector of polynomials.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> operator*(const T& a, const std::vector<P<T>>& b);

    /**
     * @brief Multiplication of a vector of polynomials with a scalar.
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] a Scalar.
     * @param[in] b Vector of polynomials.
     * @return std::vector<P<T>> Output vector of polynomials.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> operator*(const std::vector<P<T>>& b, const T& a);

    /**
     * @brief Multiplication of a vector of polynomials with a polynomial.
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] a Polynomial.
     * @param[in] b Vector of polynomials.
     * @return std::vector<P<T>> Output vector of polynomials.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> operator*(const P<T>& a, const std::vector<P<T>>& b);

    /**
     * @brief Multiplication of a vector of polynomials with a polynomial.
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] a Polynomial.
     * @param[in] b Vector of polynomials.
     * @return std::vector<P<T>> Output vector of polynomials.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> operator*(const std::vector<P<T>>& b, const P<T>& a);

    /**
     * @brief Division of a vector of polynomials with a scalar.
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] b Vector of polynomials.
     * @param[in] a Scalar.
     * @return std::vector<P<T>> Output vector.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> operator/(const std::vector<P<T>>& b, const T& a);

    /**
     * @brief Division of a vector of polynomials with a polynomial.
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] b Vector of polynomials.
     * @param[in] a Polynomial.
     * @return std::vector<P<T>> Output vector.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> operator/(const std::vector<P<T>>& b, const P<T>& a);

    #endif

}

#endif