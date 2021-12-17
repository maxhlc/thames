#ifndef THAMES_PROPAGATORS_GEQOE
#define THAMES_PROPAGATORS_GEQOE

#include "../types.h"

using namespace thames::types;

namespace thames::propagators::geqoe{

    /**
     * @brief State derivative for propagation using Generalised Equinoctial Orbital Elements (GEqOE).
     * 
     * @param geqoe GEqOE state.
     * @param geqoedot Time derivative of the GEqOE state.
     * @param t Current physical time.
     * @param mu Central body gravitational parameter.
     * @param U_func Perturbing potential function.
     * @param Ut_func Perturbing potential time derivative function.
     * @param F_func Total perturbing acceleration function.
     * @param P_func Non-potential perturbing acceleration function.
     */
    void derivative(const Vector6 &geqoe, Vector6 &geqoedot, const double t, const double &mu, const Potential &U_func, const PotentialDerivative &Ut_func, const Force &F_func, const Force &P_func);

    /**
     * @brief Propagate Cartesian state via Generalised Equinoctial Orbital Elements (GEqOE).
     * 
     * @param tstart Propagation start time in physical time.
     * @param tend Propagation end time in physical time.
     * @param tstep Initial timestep for propagation.
     * @param RV Initial Cartesian state.
     * @param mu Central body gravitational parameter.
     * @param U_func Perturbing potential function.
     * @param Ut_func Perturbing potential time derivative function.
     * @param F_func Total perturbing acceleration function.
     * @param P_func Non-potential perturbing acceleration function.
     * @param atol Solver absolute tolerance.
     * @param rtol Solver relative tolerance.
     * @return Vector6 Final Cartesian state.
     */
    Vector6 propagate(double tstart, double tend, double tstep, Vector6 RV, double mu, Potential U_func, PotentialDerivative Ut_func, Force F_func, Force P_func, double atol = 1e-10, double rtol = 1e-10);

}

#endif