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

#ifndef THAMES_CONVERSIONS_CARTESIAN
#define THAMES_CONVERSIONS_CARTESIAN

#include <vector>

namespace thames::conversions::cartesian {

    /////////////////
    // Polynomials //
    /////////////////

    /**
     * @brief Calculate lower and upper bounds for each state variable in a set of Cartesian state vectors.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-02-17
     * 
     * @tparam T Numeric type.
     * @param[in] states Cartesian state vectors. 
     * @param[out] lower Minimum of each state variable. 
     * @param[out] upper Maximum of each state variable.
     */
    template<class T>
    void states_to_bounds(const std::vector<std::vector<T>>& states, std::vector<T>& lower, std::vector<T>& upper);

    /**
     * @brief Calculate sample vector of a Cartesian state, given lower and upper bounds.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-02-18
     * 
     * @tparam T Numeric type.
     * @param[in] state Cartesian state vector. 
     * @param[in] lower Minimum of each state variable. 
     * @param[in] upper Maximum of each state variable.
     * @return std::vector<T> Sample vector.
     */
    template<class T>
    std::vector<T> state_to_sample(const std::vector<T>& state, const std::vector<T>& lower, const std::vector<T>& upper);

    /**
     * @brief Calculate sample vectors of a Cartesian state, given lower and upper bounds.
     * 
     * @tparam T Numeric type.
     * @param[in] states Cartesian state vectors. 
     * @param[in] lower Minimum of each state variable. 
     * @param[in] upper Maximum of each state variable.
     * @return std::vector<std::vector<T>> Sample vectors.
     */
    template<class T>
    std::vector<std::vector<T>> state_to_sample(const std::vector<std::vector<T>>& states, const std::vector<T>& lower, const std::vector<T>& upper);

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Generate vector of polynomials representing Cartesian state, and its uncertainty.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-01-31
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] RV Cartesian state.
     * @param[in] RVunc Cartesian state uncertainty.
     * @param[in] degree Polynomial degree.
     * @param[out] RVPolynomial Polynomial Cartesian state.
     */
    template<class T, template<class> class P>
    void cartesian_to_polynomial(const std::vector<T>& RV, const std::vector<T>& RVunc, int degree, std::vector<P<T>>& RVPolynomial);
    
    /**
     * @brief Generate vector of polynomials representing Cartesian state, and its uncertainty, from a set of Cartesian states.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-02-18
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] RVs Cartesian states.
     * @param[in] degree Polynomial degree.
     * @param[out] RVPolynomial Polynomial Cartesian state.
     * @param[out] lower Minimum of each state variable. 
     * @param[out] upper Maximum of each state variable.
     */
    template<class T, template<class> class P>
    void cartesian_to_polynomial(const std::vector<std::vector<T>>& RVs, int degree, std::vector<P<T>>& RVPolynomial, std::vector<T>& lower, std::vector<T>& upper);

    #endif

}

#endif