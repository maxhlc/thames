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

#ifndef THAMES_PROPAGATORS_COWELL
#define THAMES_PROPAGATORS_COWELL

#include <array>
#include <memory>
#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Dynamics/base_dynamics.h"
#endif

#include "basepropagator.h"
#include "../perturbations/baseperturbation.h"
#include "../constants/statetypes.h"
#include "../conversions/dimensional.h"

namespace thames::propagators {

    using thames::propagators::basepropagator::BasePropagator;
    using thames::perturbations::baseperturbation::BasePerturbation;
    using thames::conversions::dimensional::DimensionalFactors;

    ///////////
    // Reals //
    ///////////

    /**
     * @brief Propagator object for Cowell's method.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-27
     * 
     * @tparam T Numeric type.
     */
    template<class T>
    class CowellPropagator : public BasePropagator<T> {

        private:

            /// Perturbation object
            using BasePropagator<T>::m_perturbation;

            /// Dimensional factors
            using BasePropagator<T>::m_factors;

            /// Gravitational parameter
            using BasePropagator<T>::m_mu;

            /// Flag for whether to propagate in non-dimensional form
            using BasePropagator<T>::m_isNonDimensional;

            /// State type for propagation
            using BasePropagator<T>::m_propstatetype;

        public:

            /**
             * @brief Construct a new Cowell Propagator object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-27
             * 
             * @param[in] mu Gravitational parameter. 
             * @param[in] perturbation Perturbation object.
             * @param[in] factors Dimensional factors.
             */
            CowellPropagator(const T& mu, std::shared_ptr<BasePerturbation<T>> perturbation, const std::shared_ptr<DimensionalFactors<T>> factors);

            /**
             * @brief State derivative for Cowell's method propagation.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-13
             * 
             * @param[in] RV Cartesian state.
             * @param[out] RVdot Time derivative of the Cartesian state.
             * @param[in] t Current physical time.
             */
            void derivative(const std::vector<T>& RV, std::vector<T>& RVdot, const T t) const override;

    };

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using thames::propagators::basepropagator::BasePropagatorPolynomial;
    using thames::propagators::basepropagator::BasePropagatorPolynomialDynamics;
    using thames::perturbations::baseperturbation::BasePerturbationPolynomial;

    /**
     * @brief Object for Cowell's method dynamics with polynomials, compatible with the SMART-UQ schema.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-27
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     */
    template<class T, template<class> class P>
    class CowellPropagatorPolynomialDynamics : public BasePropagatorPolynomialDynamics<T, P> {

        public:

            /// Flag for whether to propagate in non-dimensional form
            using BasePropagatorPolynomialDynamics<T, P>::m_isNonDimensional;

        private:

            /// Gravitational parameter
            using BasePropagatorPolynomialDynamics<T, P>::m_mu;

            /// Perturbation object
            using BasePropagatorPolynomialDynamics<T, P>::m_perturbation;

            /// Dimensional factors
            using BasePropagatorPolynomialDynamics<T, P>::m_factors;

        public:

            /**
             * @brief Construct a new Cowell Propagator Polynomial Dynamics object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-27
             * 
             * @param[in] mu Gravitational parameter.
             * @param[in] perturbation Perturbation object.
             * @param[in] factors Dimensional factors.
             */
            CowellPropagatorPolynomialDynamics(const T& mu, const std::shared_ptr<BasePerturbationPolynomial<T, P>> perturbation, const std::shared_ptr<const DimensionalFactors<T>> factors);

            /**
             * @brief Destroy the Cowell Propagator Polynomial Dynamics object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-27
             * 
             */
            ~CowellPropagatorPolynomialDynamics();

            /**
             * @brief Evaluate the derivative of the Cowell's method dynamics.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-13
             * 
             * @param[in] t Current physical time.
             * @param[in] RV Cartesian state.
             * @param[out] RVdot Derivative of the Cartesian state.
             * @return int 
             */
            int evaluate(const T& t, const std::vector<P<T>>& RV, std::vector<P<T>>& RVdot) const override;

    };

    /**
     * @brief Propagator object for Cowell's method with polynomials.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-27
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     */
    template<class T, template<class> class P>
    class CowellPropagatorPolynomial : public BasePropagatorPolynomial<T, P> {
        
        private:

            /// Perturbation object
            using BasePropagatorPolynomial<T, P>::m_perturbation;

            /// Dimensional factors
            using BasePropagatorPolynomial<T, P>::m_factors;

            /// Gravitational parameter
            using BasePropagatorPolynomial<T, P>::m_mu;

            /// Dynamics object
            using BasePropagatorPolynomial<T, P>::m_dyn;

        public:

            /**
             * @brief Construct a new Cowell Propagator Polynomial object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-31
             * 
             * @param[in] mu Gravitational parameter.
             * @param[in] perturbation Perturbation object.
             * @param[in] factors Dimensional factors.
             */
            CowellPropagatorPolynomial(const T& mu, const std::shared_ptr<BasePerturbationPolynomial<T, P>> perturbation, const std::shared_ptr<DimensionalFactors<T>> factors);

            /**
             * @brief Destroy the Cowell Propagator Polynomial object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-27
             * 
             */
            ~CowellPropagatorPolynomial();

    };

    #endif

}

#endif