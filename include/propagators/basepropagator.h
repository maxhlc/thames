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

#ifndef THAMES_PROPAGATORS_BASEPROPAGATOR
#define THAMES_PROPAGATORS_BASEPROPAGATOR

#include <array>
#include <vector>

#include "options.h"
#include "../constants/statetypes.h"
#include "../conversions/dimensional.h"
#include "../perturbations/baseperturbation.h"

namespace thames::propagators::basepropagator {

    using thames::conversions::dimensional::DimensionalFactors;
    using thames::perturbations::baseperturbation::BasePerturbation;
    using thames::perturbations::baseperturbation::BasePerturbationPolynomial;

    ///////////
    // Reals //
    ///////////

    /**
     * @brief Base propagator abstract object.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-13
     * 
     * @tparam T Numeric type.
     */
    template<class T>
    class BasePropagator {

        protected:

            /// Perturbation object
            BasePerturbation<T>* const m_perturbation;

            /// Dimensional factors
            const DimensionalFactors<T>* m_factors;

            /// Flag for whether to propagate in non-dimensional form
            bool m_isNonDimensional = false;

        public:

            /**
             * @brief Construct a new Base Propagator object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-11
             * 
             * @param[in] perturbation Perturbation object.
             * @param[in] factors Dimensional factors.
             */
            BasePropagator(BasePerturbation<T>* const perturbation, const DimensionalFactors<T>* factors) : m_perturbation(perturbation), m_factors(factors) {

            }

            ////////////
            // Arrays //
            ////////////

            /**
             * @brief State derivative method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-26
             * 
             * @param[in] x State.
             * @param[out] dxdt State derivative.
             * @param[in] t Current time.
             */
            virtual void derivative(const std::array<T, 6>& x, std::array<T, 6>& dxdt, const T t) const;

            /**
             * @brief Propagation method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-13
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::array<T, 6> Final state.
             */
            virtual std::array<T, 6> propagate(T tstart, T tend, T tstep, std::array<T, 6> state, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype);

            /////////////
            // Vectors //
            /////////////

            /**
             * @brief State derivative method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-26
             * 
             * @param[in] x State.
             * @param[out] dxdt State derivative.
             * @param[in] t Time.
             */
            virtual void derivative(const std::vector<T>& x, std::vector<T>& dxdt, const T t) const;

            /**
             * @brief Propagation method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-13
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::vector<T> Final state.
             */
            virtual std::vector<T> propagate(T tstart, T tend, T tstep, std::vector<T> state, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype);

            /**
             * @brief Propagation method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-13
             * 
             * @param[in] tvector Vector of times.
             * @param[in] tstep Initial timestep for propagation
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::vector<std::vector<T>> Propagated states.
             */
            virtual std::vector<std::vector<T>> propagate(std::vector<T> tvector, T tstep, std::vector<T> state, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype);         

    };

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Base propagator abstract object for polynomial propagations.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-12
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     */
    template<class T, template<class> class P>
    class BasePropagatorPolynomial {

        protected:

            /// Perturbation object
            BasePerturbationPolynomial<T, P>* const m_perturbation;

            /// Dimensional factors
            const DimensionalFactors<T>* m_factors;

        public:

            /**
             * @brief Construct a new Base Propagator Polynomial object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-12
             * 
             * @param[in] perturbation Perturbation object.
             * @param[in] factors Dimensional factors.
             */
            BasePropagatorPolynomial(BasePerturbationPolynomial<T, P>* const perturbation, const DimensionalFactors<T>* factors) : m_perturbation(perturbation), m_factors(factors) {

            }

            /**
             * @brief Propagation method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-13
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::vector<P<T>> Final state.
             */
            virtual std::vector<P<T>> propagate(T tstart, T tend, T tstep, std::vector<P<T>> state, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype);

            /**
             * @brief Propagation method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-13
             * 
             * @param[in] tvector Vector of times.
             * @param[in] tstep Initial timestep for propagation
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::vector<std::vector<P<T>>> Propagated states.
             */
            virtual std::vector<std::vector<P<T>>> propagate(std::vector<T> tvector, T tstep, std::vector<P<T>> state, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype);

    };

    #endif

}

#endif