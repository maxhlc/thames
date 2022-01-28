#ifndef THAMES_PROPAGATORS_COWELL
#define THAMES_PROPAGATORS_COWELL

#include <array>
#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Dynamics/base_dynamics.h"
#endif

#include "basepropagator.h"
#include "../perturbations/baseperturbation.h"

using namespace thames::propagators::basepropagator;
using namespace thames::perturbations::baseperturbation;

namespace thames::propagators {

    ///////////
    // Reals //
    ///////////

    /**
     * @brief Propagator object for Cowell's method.
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
             * @param[in] RV Cartesian state.
             * @param[out] RVdot Time derivative of the Cartesian state.
             * @param[in] t Current physical time.
             */
            void derivative(const std::array<T, 6>& RV, std::array<T, 6>& RVdot, const T t) const override;

            /**
             * @brief Propagate Cartesian state using Cowell's method.
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] RV Initial Cartesian state.
             * @param[in] atol Solver absolute tolerance.
             * @param[in] rtol Solver relative tolerance.
             * @return std::array<T, 6> Final Cartesian state.
             */
            std::array<T, 6> propagate(T tstart, T tend, T tstep, std::array<T, 6> RV, T atol = 1e-10, T rtol = 1e-10) const override;

            /////////////
            // Vectors //
            /////////////

            /**
             * @brief State derivative for Cowell's method propagation.
             * 
             * @param[in] RV Cartesian state.
             * @param[out] RVdot Time derivative of the Cartesian state.
             * @param[in] t Current physical time.
             */
            void derivative(const std::vector<T>& RV, std::vector<T>& RVdot, const T t) const override;

            /**
             * @brief Propagate Cartesian state using Cowell's method.
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] RV Initial Cartesian state.
             * @param[in] atol Solver absolute tolerance.
             * @param[in] rtol Solver relative tolerance.
             * @return std::vector<T> Final Cartesian state.
             */
            std::vector<T> propagate(T tstart, T tend, T tstep, std::vector<T> RV, T atol = 1e-10, T rtol = 1e-10) const override;

    };

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Object for Cowell's method dynamics with polynomials, compatible with the SMART-UQ schema.
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
             * @param[in] mu Gravitational parameter.
             * @param[in] perturbation Perturbation object.
             */
            CowellPropagatorPolynomialDynamics(const T& mu, const BasePerturbationPolynomial<T, P>* perturbation);

            /**
             * @brief Destroy the Cowell Propagator Polynomial Dynamics object.
             * 
             */
            ~CowellPropagatorPolynomialDynamics();

            /**
             * @brief Evaluate the derivative of the Cowell's method dynamics.
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
             * @param[in] mu Gravitational parameter.
             * @param[in] perturbation Perturbation object.
             */
            CowellPropagatorPolynomial(const T& mu, const BasePerturbationPolynomial<T, P>* perturbation);

            /**
             * @brief Destroy the Cowell Propagator Polynomial object
             * 
             */
            ~CowellPropagatorPolynomial();

            /**
             * @brief Propagate Cartesian state using Cowell's method.
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] RV Initial Cartesian state.
             * @return std::vector<P<T>> Final Cartesian state.
             */
            std::vector<P<T>> propagate(T tstart, T tend, T tstep, std::vector<P<T>> RV) const override;

    };

    #endif

}

#endif