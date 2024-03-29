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
#include <memory>
#include <vector>
#include <string>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Dynamics/base_dynamics.h"
#endif

#include "../constants/statetypes.h"
#include "../conversions/dimensional.h"
#include "../perturbations/baseperturbation.h"
#include "../settings/settings.h"

namespace thames::propagators::basepropagator {

    using thames::constants::statetypes::StateTypes;
    using thames::conversions::dimensional::DimensionalFactors;
    using thames::perturbations::baseperturbation::BasePerturbation;
    using thames::settings::PropagatorParameters;

    ///////////
    // Reals //
    ///////////

    /**
     * @brief Base propagator object.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     */
    template<class T>
    class BasePropagator {

        protected:

            /// Gravitational parameter
            const T m_mu;

            /// Perturbation object
            const std::shared_ptr<BasePerturbation<T>> m_perturbation;

            /// Dimensional factors
            const std::shared_ptr<DimensionalFactors<T>> m_factors;

            /// Flag for whether to propagate in non-dimensional form
            bool m_isNonDimensional = false;

            /// State type for propagation
            const StateTypes m_propstatetype;

        public:

            /**
             * @brief Construct a new Base Propagator object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-27
             * 
             * @param[in] mu Gravitational parameter.
             * @param[in] perturbation Perturbation object.
             * @param[in] factors Dimensional factors.
             * @param[in] propstatetype Type of state used during propagation.
             */
            BasePropagator(const T& mu, const std::shared_ptr<BasePerturbation<T>> perturbation, const std::shared_ptr<DimensionalFactors<T>> factors, const StateTypes propstatetype);

            /**
             * @brief Destroy the Base Propagator object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-16
             * 
             */
            ~BasePropagator();

            /**
             * @brief State derivative method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-13
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
             * @date 2022-11-04
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::vector<T> Final state.
             */
            std::vector<T> propagate(T tstart, T tend, T tstep, std::vector<T> state, const PropagatorParameters<T> options, const StateTypes statetype);

            /**
             * @brief Propagation method (with intermediate output).
             * 
             * @note A separate propagation is called for each intermediate output interval, therefore any required state conversions occur multiple times.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-07-06
             * 
             * @param[in] tvec Vector of physical propagation times.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::vector<std::vector<T>> Final state.
             */
            std::vector<std::vector<T>> propagate(const std::vector<T> tvec, const T tstep, const std::vector<T> state, const PropagatorParameters<T> options, const StateTypes statetype);

            /**
             * @brief Propagation method for sets.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-06-02
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::vector<std::vector<T>> Final state.
             */
            std::vector<std::vector<T>> propagate(const T tstart, const T tend, const T tstep, const std::vector<std::vector<T>> state, const PropagatorParameters<T> options, const StateTypes statetype);

            /**
             * @brief Propagation method for sets (with intermediate output).
             * 
             * @note A separate propagation is called for each intermediate output interval, therefore any required state conversions occur multiple times.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-07-06
             * 
             * @param[in] tvec Vector of physical propagation times.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::vector<std::vector<std::vector<T>>> Final state.
             */
            std::vector<std::vector<std::vector<T>>> propagate(const std::vector<T> tvec, const T tstep, const std::vector<std::vector<T>> state, const PropagatorParameters<T> options, const StateTypes statetype);

    };

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using thames::perturbations::baseperturbation::BasePerturbationPolynomial;

    /**
     * @brief Base object for dynamics with polynomials, compatible with the SMART-UQ schema.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-27
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     */
    template<class T, template<class> class P>
    class BasePropagatorPolynomialDynamics : public smartuq::dynamics::base_dynamics<P<T>> {

        public:

            /// Flag for whether to propagate in non-dimensional form
            bool m_isNonDimensional = false;

        protected:

            /// Dynamics name
            using smartuq::dynamics::base_dynamics<P<T>>::m_name;

            /// Gravitational parameter
            const T m_mu;

            /// Perturbation object
            const std::shared_ptr<BasePerturbationPolynomial<T, P>> m_perturbation;

            /// Dimensional factors
            const std::shared_ptr<const DimensionalFactors<T>> m_factors;

        public:

            /**
             * @brief Construct a new Base Propagator Polynomial Dynamics object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-27
             * 
             * @param[in] name Dynamics object name.
             * @param[in] mu Gravitational parameter.
             * @param[in] perturbation Perturbation object.
             * @param[in] factors Dimensional factors.
             */
            BasePropagatorPolynomialDynamics(std::string name, const T& mu, const std::shared_ptr<BasePerturbationPolynomial<T, P>> perturbation, const std::shared_ptr<const DimensionalFactors<T>> factors);

            /**
             * @brief Destroy the Base Propagator Polynomial Dynamics object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-13
             * 
             */
            ~BasePropagatorPolynomialDynamics();

            /**
             * @brief Evaluate the derivative of the dynamics.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-16
             * 
             * @param[in] t Current physical time.
             * @param[in] x Current state.
             * @param[out] dxdt Derivative of the state.
             * @return int 
             */
            virtual int evaluate(const T& t, const std::vector<P<T>>& x, std::vector<P<T>>& dxdt) const;

    };

    /**
     * @brief Base propagator abstract object for polynomial propagations.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-09-28
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     */
    template<class T, template<class> class P>
    class BasePropagatorPolynomial {

        protected:

            /// Gravitational parameter
            const T m_mu;

            /// Perturbation object
            const std::shared_ptr<BasePerturbationPolynomial<T, P>> m_perturbation;

            /// Dimensional factors
            const std::shared_ptr<DimensionalFactors<T>> m_factors;

            /// Dynamics object
            const std::shared_ptr<BasePropagatorPolynomialDynamics<T, P>> m_dyn;

            /// State type for propagation
            const StateTypes m_propstatetype;

        public:

            /**
             * @brief Construct a new Base Propagator Polynomial object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-31
             * 
             * @param[in] mu Gravitational parameter.
             * @param[in] perturbation Perturbation object.
             * @param[in] factors Dimensional factors.
             * @param[in] dyn Dynamics object.
             * @param[in] propstatetype Type of state used during propagation.
             */
            BasePropagatorPolynomial(const T& mu, const std::shared_ptr<BasePerturbationPolynomial<T, P>> perturbation, const std::shared_ptr<DimensionalFactors<T>> factors, const std::shared_ptr<BasePropagatorPolynomialDynamics<T, P>> dyn, const StateTypes propstatetype);

            /**
             * @brief Propagation method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-11-04
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::vector<P<T>> Final state.
             */
            std::vector<P<T>> propagate(T tstart, T tend, T tstep, std::vector<P<T>> state, const PropagatorParameters<T> options, const StateTypes statetype);

            /**
             * @brief Propagation method (with intermediate output).
             * 
             * @note A separate propagation is called for each intermediate output interval, therefore any required state conversions occur multiple times.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-11-04
             * 
             * @param[in] tvec Vector of physical propagation times.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::vector<std::vector<P<T>>> Final state.
             */
            std::vector<std::vector<P<T>>> propagate(const std::vector<T> tvec, const T tstep, const std::vector<P<T>> state, const PropagatorParameters<T> options, const StateTypes statetype);

            /**
             * @brief Propagation method for sets of points.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-09-28
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] states Initial states.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @param[in] degree Polynomial degree.
             * @return std::vector<std::vector<T>> Final state.
             */
            std::vector<std::vector<T>> propagate(const T tstart, const T tend, const T tstep, std::vector<std::vector<T>> states, const PropagatorParameters<T> options, const StateTypes statetype, const unsigned int degree);

            /**
             * @brief Propagation method for sets of points (with intermediate output).
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-11-04
             * 
             * @param[in] tvec Vector of physical propagation times.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] states Initial states.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @param[in] degree Polynomial degree.
             * @return std::vector<std::vector<std::vector<T>>> Final state.
             */
            std::vector<std::vector<std::vector<T>>> propagate(const std::vector<T> tvec, const T tstep, const std::vector<std::vector<T>> states, const PropagatorParameters<T> options, const StateTypes statetype, const unsigned int degree);

    };

    #endif

}

#endif