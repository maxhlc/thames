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

#ifndef THAMES_CONVERSIONS_GEQOE
#define THAMES_CONVERSIONS_GEQOE

#include <array>
#include <memory>
#include <vector>

#include "../perturbations/baseperturbation.h"

namespace thames::conversions::geqoe{

    using namespace thames::perturbations::baseperturbation;

    ////////////
    // Arrays //
    ////////////

    /**
     * @brief Convert from Cartesian state to Generalised Equinoctial Orbital Elements (GEqOE).
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-26
     * 
     * @tparam T Numeric type.
     * @param[in] t Current physical time.
     * @param[in] RV Cartesian state.
     * @param[in] mu Gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @return std::array<T, 6> GEqOE state.
     */
    template<class T>
    std::array<T, 6> cartesian_to_geqoe(const T& t, const std::array<T, 6>& RV, const T& mu, const std::shared_ptr<const BasePerturbation<T>> perturbation);

    /**
     * @brief Convert from Generalised Equinoctial Orbital Elements (GEqOE) to Cartesian state.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-26
     * 
     * @tparam T Numeric type.
     * @param[in] t Current physical time.
     * @param[in] geqoe GEqOE state.
     * @param[in] mu Gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @return std::array<T, 6> Cartesian state.
     */
    template<class T>
    std::array<T, 6> geqoe_to_cartesian(const T& t, const std::array<T, 6>& geqoe, const T& mu, const std::shared_ptr<const BasePerturbation<T>> perturbation);

    /////////////
    // Vectors //
    /////////////

    /**
     * @brief Convert from Cartesian state to Generalised Equinoctial Orbital Elements (GEqOE).
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-26
     * 
     * @tparam T Numeric type.
     * @param[in] t Current physical time.
     * @param[in] RV Cartesian state.
     * @param[in] mu Gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @return std::vector<T> GEqOE state.
     */
    template<class T>
    std::vector<T> cartesian_to_geqoe(const T& t, const std::vector<T>& RV, const T& mu, const std::shared_ptr<const BasePerturbation<T>> perturbation);

    /**
     * @brief Convert from Generalised Equinoctial Orbital Elements (GEqOE) to Cartesian state.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-26
     * 
     * @tparam T Numeric type.
     * @param[in] t Current physical time.
     * @param[in] geqoe GEqOE state.
     * @param[in] mu Gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @return std::vector<T> Cartesian state.
     */
    template<class T>
    std::vector<T> geqoe_to_cartesian(const T& t, const std::vector<T>& geqoe, const T& mu, const std::shared_ptr<const BasePerturbation<T>> perturbation);

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Convert from Cartesian state to Generalised Equinoctial Orbital Elements (GEqOE).
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-26
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] t Current physical time.
     * @param[in] geqoe GEqOE state.
     * @param[in] mu Gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @return std::vector<P<T>> GEqOE state.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> cartesian_to_geqoe(const T& t, const std::vector<P<T>>& RV, const T& mu, const std::shared_ptr<const BasePerturbationPolynomial<T, P>> perturbation);

    /**
     * @brief Convert from Generalised Equinoctial Orbital Elements (GEqOE) to Cartesian state.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-26
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] t Current physical time.
     * @param[in] geqoe GEqOE state.
     * @param[in] mu Gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @return std::vector<P<T>> Cartesian state.
     */
    template<class T, template<class> class P>
    std::vector<P<T>> geqoe_to_cartesian(const T& t, const std::vector<P<T>>& geqoe, const T& mu, const std::shared_ptr<const BasePerturbationPolynomial<T, P>> perturbation);

    #endif

}

#endif