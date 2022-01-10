#ifndef THAMES_PROPAGATORS_COWELL
#define THAMES_PROPAGATORS_COWELL

#include "../types.h"

using namespace thames::types;

namespace thames::propagators::cowell{

    /**
     * @brief State derivative for Cowell's method propagation.
     * 
     * @param[in] RV Cartesian state.
     * @param[out] RVdot Time derivative of the Cartesian state.
     * @param[in] t Current physical time.
     * @param[in] mu Central body gravitational parameter.
     * @param[in] F_func Perturbing acceleration function.
     */
    void derivative(const Vector6 &RV, Vector6 &RVdot, const double t, const double &mu, const AccelerationFunc &F_func);

    /**
     * @brief Propagate Cartesian state using Cowell's method.
     * 
     * @param[in] tstart Propagation start time in physical time.
     * @param[in] tend Propagation end time in physical time.
     * @param[in] tstep Initial timestep for propagation.
     * @param[in] RV Initial Cartesian state.
     * @param[in] mu Central body gravitational parameter.
     * @param[in] F_func Perturbing acceleration function.
     * @param[in] atol Solver absolute tolerance.
     * @param[in] rtol Solver relative tolerance.
     * @return Vector6 Final Cartesian state.
     */
    Vector6 propagate(double tstart, double tend, double tstep, Vector6 RV, double mu, AccelerationFunc F_func, double atol = 1e-10, double rtol = 1e-10);

}

#endif