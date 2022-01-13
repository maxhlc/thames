#ifndef THAMES_PROPAGATORS_GEQOE
#define THAMES_PROPAGATORS_GEQOE

#include "../perturbations/baseperturbation.h"
#include "../types.h"

using namespace thames::types;
using namespace thames::perturbations::baseperturbation;

namespace thames::propagators::geqoe{

    /**
     * @brief State derivative for propagation using Generalised Equinoctial Orbital Elements (GEqOE).
     * 
     * @tparam real Type for real numbers.
     * @tparam vector3 Type for the state vector slices.
     * @tparam vector6 Type for the state vector.
     * @param[in] geqoe GEqOE state.
     * @param[in,out] geqoedot Time derivative of the GEqOE state.
     * @param[in] t Current physical time.
     * @param[in] mu Central body gravitational parameter.
     * @param[in] perturbation Perturbation object.
     */
    template<class real, class vector3, class vector6>
    void derivative(const vector6 &geqoe, vector6 &geqoedot, const real t, const real &mu, BasePerturbation<real, vector3> &perturbation);

    /**
     * @brief Propagate Cartesian state via Generalised Equinoctial Orbital Elements (GEqOE).
     * 
     * @tparam real Type for real numbers.
     * @tparam vector3 Type for the state vector slices.
     * @tparam vector6 Type for the state vector.
     * @param[in] tstart Propagation start time in physical time.
     * @param[in] tend Propagation end time in physical time.
     * @param[in] tstep Initial timestep for propagation.
     * @param[in] RV Initial Cartesian state.
     * @param[in] mu Central body gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @param[in] atol Solver absolute tolerance.
     * @param[in] rtol Solver relative tolerance.
     * @return Vector6 Final Cartesian state.
     */
    template<class real, class vector3, class vector6>
    vector6 propagate(real tstart, real tend, real tstep, vector6 RV, real mu, BasePerturbation<real, vector3> &perturbation, real atol = 1e-10, real rtol = 1e-10);

}

#endif