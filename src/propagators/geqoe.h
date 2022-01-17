#ifndef THAMES_PROPAGATORS_GEQOE
#define THAMES_PROPAGATORS_GEQOE

#include <array>

#include "../perturbations/baseperturbation.h"

using namespace thames::perturbations::baseperturbation;

namespace thames::propagators::geqoe{

    /**
     * @brief State derivative for propagation using Generalised Equinoctial Orbital Elements (GEqOE).
     * 
     * @tparam T Numeric type.
     * @param[in] geqoe GEqOE state.
     * @param[in,out] geqoedot Time derivative of the GEqOE state.
     * @param[in] t Current physical time.
     * @param[in] mu Central body gravitational parameter.
     * @param[in] perturbation Perturbation object.
     */
    template<class T>
    void derivative(const std::array<T, 6>& geqoe, std::array<T, 6>& geqoedot, const T t, const T& mu, const BasePerturbation<T>& perturbation);

    /**
     * @brief Propagate Cartesian state via Generalised Equinoctial Orbital Elements (GEqOE).
     * 
     * @tparam T Numeric type.
     * @param[in] tstart Propagation start time in physical time.
     * @param[in] tend Propagation end time in physical time.
     * @param[in] tstep Initial timestep for propagation.
     * @param[in] RV Initial Cartesian state.
     * @param[in] mu Central body gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @param[in] atol Solver absolute tolerance.
     * @param[in] rtol Solver relative tolerance.
     * @return std::array<T, 6> Final Cartesian state.
     */
    template<class T>
    std::array<T, 6> propagate(T tstart, T tend, T tstep, std::array<T, 6> RV, T mu, const BasePerturbation<T>& perturbation, T atol = 1e-10, T rtol = 1e-10);

}

#endif