#ifndef THAMES_PROPAGATORS_GEQOE
#define THAMES_PROPAGATORS_GEQOE

#include <array>
#include <vector>

#include "basepropagator.h"
#include "../perturbations/baseperturbation.h"

using namespace thames::propagators::basepropagator;
using namespace thames::perturbations::baseperturbation;

namespace thames::propagators {

    template<class T>
    class GEqOEPropagator : public BasePropagator<T> {

        private:

            /// Gravitational parameter
            const T m_mu;

            /// Perturbation object
            const BasePerturbation<T>* m_perturbation;

        public:

            /**
             * @brief Construct a new GEqOE Propagator object.
             * 
             * @param[in] mu Gravitational parameter.
             * @param[in] perturbation Perturbation object. 
             */
            GEqOEPropagator(const T& mu, const BasePerturbation<T>* perturbation);

            ////////////
            // Arrays //
            ////////////

            /**
             * @brief State derivative for propagation using Generalised Equinoctial Orbital Elements (GEqOE).
             * 
             * @param[in] geqoe GEqOE state.
             * @param[out] geqoedot Time derivative of the GEqOE state.
             * @param[in] t Current physical time.
             */
            void derivative(const std::array<T, 6>& geqoe, std::array<T, 6>& geqoedot, const T t) const override;

            /**
             * @brief Propagate Cartesian state via Generalised Equinoctial Orbital Elements (GEqOE).
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
             * @brief State derivative for propagation using Generalised Equinoctial Orbital Elements (GEqOE).
             * 
             * @param[in] geqoe GEqOE state.
             * @param[out] geqoedot Time derivative of the GEqOE state.
             * @param[in] t Current physical time.
             */
            void derivative(const std::vector<T>& geqoe, std::vector<T>& geqoedot, const T t) const override;

            /**
             * @brief Propagate Cartesian state via Generalised Equinoctial Orbital Elements (GEqOE).
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