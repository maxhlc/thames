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
     * @param[in] geqoe GEqOE state.
     * @param[in,out] geqoedot Time derivative of the GEqOE state.
     * @param[in] t Current physical time.
     * @param[in] mu Central body gravitational parameter.
     * @param[in] perturbation Perturbation object.
     */
    void derivative(const Vector6 &geqoe, Vector6 &geqoedot, const double t, const double &mu, BasePerturbation<double, Vector3> &perturbation);

    /**
     * @brief Propagate Cartesian state via Generalised Equinoctial Orbital Elements (GEqOE).
     * 
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
    Vector6 propagate(double tstart, double tend, double tstep, Vector6 RV, double mu, BasePerturbation<double, Vector3> &perturbation, double atol = 1e-10, double rtol = 1e-10);

}

#endif