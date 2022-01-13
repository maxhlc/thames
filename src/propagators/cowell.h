#ifndef THAMES_PROPAGATORS_COWELL
#define THAMES_PROPAGATORS_COWELL

#include "../perturbations/baseperturbation.h"
#include "../types.h"

using namespace thames::types;
using namespace thames::perturbations::baseperturbation;

namespace thames::propagators::cowell{

    /**
     * @brief State derivative for Cowell's method propagation.
     * 
     * @tparam real Type for real numbers.
     * @tparam vector3 Type for the state vector slices.
     * @tparam vector6 Type for the state vector.
     * @param[in] RV Cartesian state.
     * @param[out] RVdot Time derivative of the Cartesian state.
     * @param[in] t Current physical time.
     * @param[in] mu Central body gravitational parameter.
     * @param[in] perturbation Perturbation object.
     */
    template<class real, class vector3, class vector6>
    void derivative(const vector6 &RV, vector6 &RVdot, const real t, const real &mu, BasePerturbation<real, vector3> &perturbation);

    /**
     * @brief Propagate Cartesian state using Cowell's method.
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