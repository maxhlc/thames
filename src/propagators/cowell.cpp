#include <cmath>
#include <functional>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

#include "cowell.h"
#include "../perturbations/baseperturbation.h"
#include "../util/vector.h"
#include "../types.h"

using namespace thames::types;
using namespace thames::perturbations::baseperturbation;

namespace thames::propagators::cowell{

    template<class real, class vector3, class vector6>
    void derivative(const vector6 &RV, vector6 &RVdot, const real t, const real &mu, BasePerturbation<real, vector3> &perturbation) {
        // Extract Cartesian state vectors
        vector3 R, V;
        R = thames::util::vector::slice<vector3, vector6, unsigned int>(RV, 0, 2);
        V = thames::util::vector::slice<vector3, vector6, unsigned int>(RV, 3, 5);

        // Calculate range
        real r = thames::util::vector::norm3<real, vector3>(R);

        // Calculate perturbing acceleration
        vector3 F = perturbation.acceleration_total(t, R, V);

        // Calculate central body acceleration
        vector3 G = thames::util::vector::mult3<real, vector3>(-mu/pow(r, 3.0), R);

        // Calculate acceleration
        vector3 A;
        for(unsigned int ii=0; ii<3; ii++)
            A[ii] = G[ii] + F[ii];

        // Store state derivative
        for(unsigned int ii=0; ii<3; ii++){
            RVdot[ii] = V[ii];
            RVdot[ii+3] = A[ii];
        }
    }
    template void derivative<double, Vector3, Vector6>(const Vector6&, Vector6&, const double, const double&, BasePerturbation<double, Vector3>&);

    template<class real, class vector3, class vector6>
    vector6 propagate(real tstart, real tend, real tstep, vector6 RV, real mu, BasePerturbation<real, vector3> &perturbation, real atol, real rtol){
        // Declare derivative function wrapper
        auto derivative_param = [&](const vector6 &x, vector6 &dxdt, const real time){
            derivative<double, Vector3, Vector6>(x, dxdt, time, mu, perturbation);
        };

        // Declare stepper
        boost::numeric::odeint::runge_kutta_cash_karp54<vector6> stepper;
        auto steppercontrolled = boost::numeric::odeint::make_controlled(atol, rtol, stepper);

        // Propagate orbit
        boost::numeric::odeint::integrate_adaptive(steppercontrolled, derivative_param, RV, tstart, tend, tstep);

        // Return final state
        return RV;
    }
    template Vector6 propagate<double, Vector3, Vector6>(double, double, double, Vector6, double, BasePerturbation<double, Vector3>&, double, double);

}