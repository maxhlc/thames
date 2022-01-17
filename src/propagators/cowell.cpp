#include <array>
#include <cmath>
#include <functional>

#include <boost/numeric/odeint.hpp>

#include "cowell.h"
#include "../perturbations/baseperturbation.h"
#include "../util/vector.h"

using namespace thames::perturbations::baseperturbation;

namespace thames::propagators::cowell{

    template<class T>
    void derivative(const std::array<T, 6>& RV, std::array<T, 6>& RVdot, const T t, const T& mu, const BasePerturbation<T>& perturbation) {
        // Extract Cartesian state vectors
        std::array<T, 3> R, V;
        R = thames::util::vector::slice<T, 6, 3>(RV, 0, 2);
        V = thames::util::vector::slice<T, 6, 3>(RV, 3, 5);

        // Calculate range
        T r = thames::util::vector::norm3<T>(R);

        // Calculate perturbing acceleration
        std::array<T, 3> F = perturbation.acceleration_total(t, R, V);

        // Calculate central body acceleration
        std::array<T, 3> G = thames::util::vector::mult3<T>(-mu/pow(r, 3.0), R);

        // Calculate acceleration
        std::array<T, 3> A;
        for(unsigned int ii=0; ii<3; ii++)
            A[ii] = G[ii] + F[ii];

        // Store state derivative
        for(unsigned int ii=0; ii<3; ii++){
            RVdot[ii] = V[ii];
            RVdot[ii+3] = A[ii];
        }
    }
    template void derivative<double>(const std::array<double, 6>&, std::array<double, 6>&, const double, const double&, const BasePerturbation<double>&);

    template<class T>
    std::array<T, 6> propagate(T tstart, T tend, T tstep, std::array<T, 6> RV, T mu, const BasePerturbation<T>& perturbation, T atol, T rtol){
        // Declare derivative function wrapper
        auto derivative_param = [&](const std::array<T, 6>& x, std::array<T, 6>& dxdt, const T time){
            derivative<double>(x, dxdt, time, mu, perturbation);
        };

        // Declare stepper
        boost::numeric::odeint::runge_kutta_cash_karp54<std::array<T, 6>> stepper;
        auto steppercontrolled = boost::numeric::odeint::make_controlled(atol, rtol, stepper);

        // Propagate orbit
        boost::numeric::odeint::integrate_adaptive(steppercontrolled, derivative_param, RV, tstart, tend, tstep);

        // Return final state
        return RV;
    }
    template std::array<double, 6> propagate<double>(double, double, double, std::array<double, 6>, double, const BasePerturbation<double>&, double, double);

}