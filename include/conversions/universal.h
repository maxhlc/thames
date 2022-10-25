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
#include <memory>
#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../constants/statetypes.h"
#include "../../include/conversions/dimensional.h"
#include "../perturbations/baseperturbation.h"

namespace thames::conversions::universal {

    using thames::constants::statetypes::StateTypes;
    using thames::conversions::dimensional::DimensionalFactors;
    using thames::perturbations::baseperturbation::BasePerturbation;

    ///////////
    // Reals //
    ///////////

    /**
     * @brief Universal state conversion.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-26
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
    std::vector<T> convert_state(const T& t, const std::vector<T>& state, const T& mu, const StateTypes& statetype1, const StateTypes& statetype2, const std::shared_ptr<const BasePerturbation<T>> perturbation);

    /**
     * @brief Universal state non-dimensionalisation.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @param[in] state State.
     * @param[in] statetype State type.
     * @param[in] factors Structure containing the factors for non-dimensionalisation.
     * @return std::vector<T> Non-dimensional state.
     */
    template<class T>
    std::vector<T> nondimensionalise_state(const std::vector<T>& state, const StateTypes& statetype, const DimensionalFactors<T>& factors);

    /**
     * @brief Universal state dimensionalisation.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @param[in] statend State.
     * @param[in] statetype State type.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     * @return std::vector<T> Dimensional state.
     */
    template<class T>
    std::vector<T> dimensionalise_state(const std::vector<T>& statend, const StateTypes& statetype, const DimensionalFactors<T>& factors);

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using thames::perturbations::baseperturbation::BasePerturbationPolynomial;

    /**
     * @brief Universal state conversion.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-26
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
    std::vector<P<T>> convert_state(const T& t, const std::vector<P<T>>& state, const T& mu, const StateTypes& statetype1, const StateTypes& statetype2, const std::shared_ptr<const BasePerturbationPolynomial<T, P>> perturbation);

    /**
     * @brief Universal state non-dimensionalisation.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] state State.
     * @param[in] statetype State type.
     * @param[in] factors Structure containing the factors for non-dimensionalisation.
     * @return std::vector<P<T>> Non-dimensional state.
     */
    template<class T, template <class> class P>
    std::vector<P<T>> nondimensionalise_state(const std::vector<P<T>>& state, const StateTypes& statetype, const DimensionalFactors<T>& factors);

    /**
     * @brief Universal state dimensionalisation.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] statend Non-dimensional state.
     * @param[in] statetype State type.
     * @param[in] factors Structure containing the factors for dimensionalisation.
     * @return std::vector<P<T>> Dimensional state.
     */
    template<class T, template <class> class P>
    std::vector<P<T>> dimensionalise_state(const std::vector<P<T>>& statend, const StateTypes& statetype, const DimensionalFactors<T>& factors);

    #endif

}

#endif