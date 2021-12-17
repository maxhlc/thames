#ifndef THAMES_PROPAGATORS_COWELL
#define THAMES_PROPAGATORS_COWELL

#include "../types.h"

using namespace thames::types;

namespace thames::propagators::cowell{

    /**
     * @brief State derivative for Cowell's method propagation.
     * 
     * @param RV Cartesian state.
     * @param RVdot Time derivative of the Cartesian state.
     * @param t Current physical time.
     * @param mu Central body gravitational parameter.
     * @param F_func Perturbing acceleration function.
     */
    void derivative(const Vector6 &RV, Vector6 &RVdot, const double t, const double &mu, const Force &F_func);

    /**
     * @brief Propagate Cartesian state using Cowell's method.
     * 
     * @param tstart Propagation start time in physical time.
     * @param tend Propagation end time in physical time.
     * @param tstep Initial timestep for propagation.
     * @param RV Initial Cartesian state.
     * @param mu Central body gravitational parameter.
     * @param F_func Perturbing acceleration function.
     * @param atol Solver absolute tolerance.
     * @param rtol Solver relative tolerance.
     * @return Vector6 Final Cartesian state.
     */
    Vector6 propagate(double tstart, double tend, double tstep, Vector6 RV, double mu, Force F_func, double atol = 1e-10, double rtol = 1e-10);

}

#endif