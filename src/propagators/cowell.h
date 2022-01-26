#ifndef THAMES_PROPAGATORS_COWELL
#define THAMES_PROPAGATORS_COWELL

#include <array>
#include <vector>

#include "basepropagator.h"
#include "../perturbations/baseperturbation.h"

using namespace thames::propagators::basepropagator;
using namespace thames::perturbations::baseperturbation;

namespace thames::propagators {

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

}

#endif