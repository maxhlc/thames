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

#ifndef THAMES_CONVERSIONS_UNIVERSAL
#define THAMES_CONVERSIONS_UNIVERSAL

#include <array>
#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../constants/statetypes.h"
#include "../perturbations/baseperturbation.h"

namespace thames::conversions::universal {

    using thames::constants::statetypes::StateTypes;
    using thames::perturbations::baseperturbation::BasePerturbation;

    ////////////
    // Arrays //
    ////////////

    /**
     * @brief Universal state conversion.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-13
     * 
     * @tparam T Numeric type.
     * @param[in] t Time.
     * @param[in] state State.
     * @param[in] mu Gravitational parameter.
     * @param[in] statetype1 Input state type.
     * @param[in] statetype2 Output state type.
     * @param[in] perturbation Perturbation object.
     * @return std::array<T, 6> Output state.
     */
    template<class T>
    std::array<T, 6> convert_state(const T& t, const std::array<T, 6>& state, const T& mu, const StateTypes& statetype1, const StateTypes& statetype2, const BasePerturbation<T>* const perturbation);

    /////////////
    // Vectors //
    /////////////

    /**
     * @brief Universal state conversion.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-13
     * 
     * @tparam T Numeric type.
     * @param[in] t Time.
     * @param[in] state State.
     * @param[in] mu Gravitational parameter.
     * @param[in] statetype1 Input state type.
     * @param[in] statetype2 Output state type.
     * @param[in] perturbation Perturbation object.
     * @return std::vector<T> Output state.
     */
    template<class T>
    std::vector<T> convert_state(const T& t, const std::vector<T>& state, const T& mu, const StateTypes& statetype1, const StateTypes& statetype2, const BasePerturbation<T>* const perturbation);

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using thames::perturbations::baseperturbation::BasePerturbationPolynomial;

    /**
     * @brief Universal state conversion.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-13
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] t Time.
     * @param[in] state State.
     * @param[in] mu Gravitational parameter.
     * @param[in] statetype1 Input state type.
     * @param[in] statetype2 Output state type.
     * @param[in] perturbation Perturbation object.
     * @return std::vector<P<T>> Output state.
     */
    template<class T, template <class> class P>
    std::vector<P<T>> convert_state(const T& t, const std::vector<P<T>>& state, const T& mu, const StateTypes& statetype1, const StateTypes& statetype2, const BasePerturbationPolynomial<T, P>* const perturbation);

    #endif

}

#endif