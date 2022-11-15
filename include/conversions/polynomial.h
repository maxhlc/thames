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

#ifndef THAMES_CONVERSIONS_POLYNOMIAL
#define THAMES_CONVERSIONS_POLYNOMIAL

namespace thames::conversions::polynomial {

    /**
     * @brief Calculate lower and upper bounds for each state variable in a set of state vectors.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @param[in] states State vectors. 
     * @param[out] lower Minimum of each state variable. 
     * @param[out] upper Maximum of each state variable.
     */
    template<class T>
    void states_to_bounds(const std::vector<std::vector<T>>& states, std::vector<T>& lower, std::vector<T>& upper);

    /**
     * @brief Calculate sample vector of a state, given lower and upper bounds.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @param[in] state State vector. 
     * @param[in] lower Minimum of each state variable. 
     * @param[in] upper Maximum of each state variable.
     * @return std::vector<T> Sample vector.
     */
    template<class T>
    std::vector<T> state_to_sample(const std::vector<T>& state, const std::vector<T>& lower, const std::vector<T>& upper);

    /**
     * @brief Calculate sample vectors of states, given lower and upper bounds.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @param[in] states State vectors. 
     * @param[in] lower Minimum of each state variable. 
     * @param[in] upper Maximum of each state variable.
     * @return std::vector<std::vector<T>> Sample vectors.
     */
    template<class T>
    std::vector<std::vector<T>> state_to_sample(const std::vector<std::vector<T>>& states, const std::vector<T>& lower, const std::vector<T>& upper);

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Generate vector of polynomials representing the state, and its uncertainty.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] state State.
     * @param[in] stateunc State uncertainty.
     * @param[in] degree Polynomial degree.
     * @param[out] statepolynomial Polynomial state.
     */
    template<class T, template<class> class P>
    void states_to_polynomial(const std::vector<T>& state, const std::vector<T>& stateunc, int degree, std::vector<P<T>>& statepolynomial);
    
    /**
     * @brief Generate vector of polynomials representing the state, and its uncertainty, from a set of states.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] states States.
     * @param[in] degree Polynomial degree.
     * @param[out] statepolynomial Polynomial state.
     * @param[out] lower Minimum of each state variable. 
     * @param[out] upper Maximum of each state variable.
     */
    template<class T, template<class> class P>
    void states_to_polynomial(const std::vector<std::vector<T>>& states, int degree, std::vector<P<T>>& statepolynomial, std::vector<T>& lower, std::vector<T>& upper);

    #endif

}

#endif