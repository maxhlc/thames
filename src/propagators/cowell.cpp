#include <functional>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

#include "cowell.h"
#include "../types.h"

using namespace thames::types;

namespace thames::propagators::cowell{

    void derivative(const Vector6 &RV, Vector6 &RVdot, const double t, const double &mu, const Force &F_func) {
        // TODO: documentation

        // Extract Cartesian state vectors
        Vector3 R, V;
        R = RV(Eigen::seq(0,2));
        V = RV(Eigen::seq(3,5));

        // Calculate range
        double r = R.norm();

        // Calculate perturbing force
        Vector3 F = F_func(t, R, V);

        // Calculate central body acceleration
        Vector3 G = -mu*R/pow(r, 3.0);

        // Calculate acceleration
        Vector3 A = G + F;

        // Store state derivative
        RVdot << V, A;    
    }

    Vector6 propagate(double tstart, double tend, double tstep, Vector6 RV, double mu, Force F_func, double atol, double rtol){
        // TODO: documentation

        // Declare derivative function wrapper
        auto derivative_param = [&](const Vector6 &x, Vector6 &dxdt, const double time){
            derivative(x, dxdt, time, mu, F_func);
        };

        // Declare stepper
        boost::numeric::odeint::runge_kutta_cash_karp54<Vector6> stepper;
        auto steppercontrolled = boost::numeric::odeint::make_controlled(atol, rtol, stepper);

        // Propagate orbit
        boost::numeric::odeint::integrate_adaptive(steppercontrolled, derivative_param, RV, tstart, tend, tstep);

        // Return final state
        return RV;
    }

}