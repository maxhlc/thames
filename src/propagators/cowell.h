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
#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Dynamics/base_dynamics.h"
#endif

#include "basepropagator.h"
#include "options.h"
#include "../perturbations/baseperturbation.h"
#include "../constants/statetypes.h"

using namespace thames::propagators::basepropagator;
using namespace thames::perturbations::baseperturbation;

namespace thames::propagators {

    ///////////
    // Reals //
    ///////////

    /**
     * @brief Propagator object for Cowell's method.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-02-28
     * 
     * @tparam T Numeric type.
     */
    template<class T>
    class CowellPropagator : public BasePropagator<T> {

        private:

            /// Gravitational parameter
            const T m_mu;

            /// Perturbation object
            const BasePerturbation<T>* m_perturbation;

        public:

            /**
             * @brief Construct a new Cowell Propagator object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-25
             * 
             * @param[in] mu Gravitational parameter. 
             * @param[in] perturbation Perturbation object.
             */
            CowellPropagator(const T& mu, const BasePerturbation<T>* perturbation);

            ////////////
            // Arrays //
            ////////////

            /**
             * @brief State derivative for Cowell's method propagation.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-25
             * 
             * @param[in] RV Cartesian state.
             * @param[out] RVdot Time derivative of the Cartesian state.
             * @param[in] t Current physical time.
             */
            void derivative(const std::array<T, 6>& RV, std::array<T, 6>& RVdot, const T t) const override;

            /**
             * @brief Propagate Cartesian state using Cowell's method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-03-07
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::array<T, 6> Final state.
             */
            std::array<T, 6> propagate(T tstart, T tend, T tstep, std::array<T, 6> state, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype = thames::constants::statetypes::CARTESIAN) const override;

            /////////////
            // Vectors //
            /////////////

            /**
             * @brief State derivative for Cowell's method propagation.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-25
             * 
             * @param[in] RV Cartesian state.
             * @param[out] RVdot Time derivative of the Cartesian state.
             * @param[in] t Current physical time.
             */
            void derivative(const std::vector<T>& RV, std::vector<T>& RVdot, const T t) const override;

            /**
             * @brief Propagate Cartesian state using Cowell's method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-03-07
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::vector<T> Final state.
             */
            std::vector<T> propagate(T tstart, T tend, T tstep, std::vector<T> state, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype = thames::constants::statetypes::CARTESIAN) const override;
            
            /**
             * @brief Propagate Cartesian state using Cowell's method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-03-07
             * 
             * @param[in] tvector Vector of times.
             * @param[in] tstep Initial timestep for propagation
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::vector<std::vector<T>> Propagated states.
             */
            std::vector<std::vector<T>> propagate(std::vector<T> tvector, T tstep, std::vector<T> RV, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype = thames::constants::statetypes::CARTESIAN) const override;

    };

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Object for Cowell's method dynamics with polynomials, compatible with the SMART-UQ schema.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-01-28
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     */
    template<class T, template<class> class P>
    class CowellPropagatorPolynomialDynamics : public smartuq::dynamics::base_dynamics<P<T>> {

        private:

            /// Dynamics name
            using smartuq::dynamics::base_dynamics<P<T>>::m_name;

            /// Gravitational parameter
            const T m_mu;

            /// Perturbation object
            const BasePerturbationPolynomial<T, P>* m_perturbation;

        public:

            /**
             * @brief Construct a new Cowell Propagator Polynomial Dynamics object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-28
             * 
             * @param[in] mu Gravitational parameter.
             * @param[in] perturbation Perturbation object.
             */
            CowellPropagatorPolynomialDynamics(const T& mu, const BasePerturbationPolynomial<T, P>* perturbation);

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
             * @date 2022-01-27
             * 
             * @param[in] t Current physical time.
             * @param[in] RV Cartesian state.
             * @param[out] RVdot Derivative of the Cartesian state.
             * @return int 
             */
            int evaluate(const T& t, const std::vector<P<T>>& RV, std::vector<P<T>>& RVdot) const;

    };

    /**
     * @brief Propagator object for Cowell's method with polynomials.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-02-28
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     */
    template<class T, template<class> class P>
    class CowellPropagatorPolynomial : public BasePropagatorPolynomial<T, P> {
        
        private:

            /// Dynamics object
            const CowellPropagatorPolynomialDynamics<T, P> m_dyn;

        public:

            /**
             * @brief Construct a new Cowell Propagator Polynomial object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-28
             * 
             * @param[in] mu Gravitational parameter.
             * @param[in] perturbation Perturbation object.
             */
            CowellPropagatorPolynomial(const T& mu, const BasePerturbationPolynomial<T, P>* perturbation);

            /**
             * @brief Destroy the Cowell Propagator Polynomial object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-27
             * 
             */
            ~CowellPropagatorPolynomial();

            /**
             * @brief Propagate Cartesian state using Cowell's method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-03-07
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::vector<P<T>> Final state.
             */
            std::vector<P<T>> propagate(T tstart, T tend, T tstep, std::vector<P<T>> state, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype = thames::constants::statetypes::CARTESIAN) const override;

            /**
             * @brief Propagate Cartesian state using Cowell's method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-03-07
             * 
             * @param[in] tvector Vector of times.
             * @param[in] tstep Initial timestep for propagation
             * @param[in] state Initial state.
             * @param[in] options Propagator options.
             * @param[in] statetype State type.
             * @return std::vector<std::vector<P<T>>> Propagated states.
             */
            std::vector<std::vector<P<T>>> propagate(std::vector<T> tvector, T tstep, std::vector<P<T>> state, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype = thames::constants::statetypes::CARTESIAN) const override;

    };

    #endif

}

#endif