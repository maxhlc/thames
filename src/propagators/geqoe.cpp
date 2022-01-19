#include <array>
#include <cmath>
#include <vector>

#include <boost/numeric/odeint.hpp>

#include "geqoe.h"
#include "../conversions/geqoe.h"
#include "../perturbations/baseperturbation.h"
#include "../util/vector.h"

using namespace thames::perturbations::baseperturbation;

namespace thames::propagators::geqoe{

    template<class T>
    void derivative(const std::array<T, 6>& geqoe, std::array<T, 6>& geqoedot, const T t, const T& mu, const BasePerturbation<T>& perturbation){
        // Extract elements
        T nu = geqoe[0];
        T p1 = geqoe[1];
        T p2 = geqoe[2];
        T q1 = geqoe[4];
        T q2 = geqoe[5];

        // Calculate Cartesian state and extract vectors
        std::array<T, 6> RV = thames::conversions::geqoe::geqoe_to_cartesian<T>(t, geqoe, mu, perturbation);
        std::array<T, 3> R, V;
        R = thames::util::vector::slice<T, 6, 3>(RV, 0, 2);
        V = thames::util::vector::slice<T, 6, 3>(RV, 3, 5);

        // Calculate range and range rate
        T r = thames::util::vector::norm3<T>(R);
        T drdt = thames::util::vector::dot3<T>(R, V)/r;

        // Calculate perturbations
        T U = perturbation.potential(t, R);
        T Ut = perturbation.potential_derivative(t, R, V);
        std::array<T, 3> F = perturbation.acceleration_total(t, R, V);
        std::array<T, 3> P = perturbation.acceleration_nonpotential(t, R, V);

        // Calculate time derivative of total energy
        T edot = Ut + thames::util::vector::dot3<T>(P, V);

        // Calculate time derivative of nu
        T nudot = -3.0*pow(nu/pow(mu, 2.0), 1.0/3.0)*edot;

        // Calculate equinocital reference frame unit vectors
        T efac = 1.0/(1.0 + pow(q1, 2.0) + pow(q2, 2.0));
        std::array<T, 3> ex, ey;
        ex[0] = efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0));
        ex[1] = efac*(2.0*q1*q2);
        ex[2] = efac*(-2.0*q1);
        ey[0] = efac*(2.0*q1*q2);
        ey[1] = efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0));
        ey[2] = efac*(2.0*q2);

        // Calculate radial unit vector
        std::array<T, 3> er = thames::util::vector::mult3<T>(1.0/r, R);

        // Calculate trig of the true longitude
        T cl = thames::util::vector::dot3<T>(er, ex);
        T sl = thames::util::vector::dot3<T>(er, ey);

        // Calculate equinoctial reference frame velocity components
        T hwh = q1*cl - q2*sl;

        // Calculate angular momentum
        std::array<T, 3> H = thames::util::vector::cross3<T>(R, V);
        T h = thames::util::vector::norm3<T>(H);
        std::array<T, 3> eh = thames::util::vector::mult3<T>(1.0/h, H);

        // Calculate the effective potential energy
        T ueff = pow(h, 2.0)/(2.0*pow(r, 2.0)) + U;

        // Calculate the generalised angular momentum
        T c = sqrt(2.0*pow(r, 2.0)*ueff);

        // Calculate the generalised semi-latus rectum
        T p = pow(c, 2.0)/mu;

        // Calculate perturbation components
        T Fr = thames::util::vector::dot3<T>(F, er);
        T Fh = thames::util::vector::dot3<T>(F, eh);

        // Calculate non-dimensional quantities
        T zeta = r/p;
        T zetatilde = 1 + zeta;

        // Calculate time derivatives of the second and third elements
        T p1dot = p2*((h - c)/pow(r, 2.0) - r/h*hwh*Fh) + 1.0/c*(r*drdt/c*p1 + zetatilde*p2 + zeta*cl)*(2.0*U - r*Fr) + r/mu*(zeta*p1 + zetatilde*sl)*edot;
        T p2dot = p1*(r/h*hwh*Fh - (h-c)/pow(r, 2.0)) + 1.0/c*(r*drdt/c*p2 - zetatilde*p1 - zeta*sl)*(2.0*U - r*Fr) + r/mu*(zeta*p2 + zetatilde*cl)*edot;

        // Calculate generalised semi-major axis
        T a = pow(mu/pow(nu, 2.0), 1.0/3.0);

        // Calculate alpha
        T alpha = 1.0/(1.0 + sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0)));

        // Calculate time derivative of the generalised mean longitude
        T Ldot = nu + (h - c)/pow(r, 2.0) - r/h*hwh*Fh + (r*drdt*c/pow(mu, 2.0)*zetatilde*alpha)*edot + 1.0/c*(1.0/alpha + alpha*(1.0 - r/a))*(2.0*U - r*Fr);

        // Calculate time derivatives of the remaining elements
        T q1dot = r/(2.0*h)*Fh*(1.0 + pow(q1, 2.0) + pow(q2, 2.0))*sl;
        T q2dot = r/(2.0*h)*Fh*(1.0 + pow(q1, 2.0) + pow(q2, 2.0))*cl;

        // Store derivatives
        geqoedot[0] = nudot;
        geqoedot[1] = p1dot;
        geqoedot[2] = p2dot;
        geqoedot[3] = Ldot;
        geqoedot[4] = q1dot;
        geqoedot[5] = q2dot;
    }
    template void derivative<double>(const std::array<double, 6>&, std::array<double, 6>&, const double, const double&, const BasePerturbation<double>&);

    template<class T>
    void derivative(const std::vector<T>& geqoe, std::vector<T>& geqoedot, const T t, const T& mu, const BasePerturbation<T>& perturbation){
        // Extract elements
        T nu = geqoe[0];
        T p1 = geqoe[1];
        T p2 = geqoe[2];
        T q1 = geqoe[4];
        T q2 = geqoe[5];

        // Calculate Cartesian state and extract vectors
        std::vector<T> RV = thames::conversions::geqoe::geqoe_to_cartesian<T>(t, geqoe, mu, perturbation);
        std::vector<T> R(3), V(3);
        R = thames::util::vector::slice<T>(RV, 0, 2);
        V = thames::util::vector::slice<T>(RV, 3, 5);

        // Calculate range and range rate
        T r = thames::util::vector::norm3<T>(R);
        T drdt = thames::util::vector::dot3<T>(R, V)/r;

        // Calculate perturbations
        T U = perturbation.potential(t, R);
        T Ut = perturbation.potential_derivative(t, R, V);
        std::vector<T> F = perturbation.acceleration_total(t, R, V);
        std::vector<T> P = perturbation.acceleration_nonpotential(t, R, V);

        // Calculate time derivative of total energy
        T edot = Ut + thames::util::vector::dot3<T>(P, V);

        // Calculate time derivative of nu
        T nudot = -3.0*pow(nu/pow(mu, 2.0), 1.0/3.0)*edot;

        // Calculate equinocital reference frame unit vectors
        T efac = 1.0/(1.0 + pow(q1, 2.0) + pow(q2, 2.0));
        std::vector<T> ex(3), ey(3);
        ex[0] = efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0));
        ex[1] = efac*(2.0*q1*q2);
        ex[2] = efac*(-2.0*q1);
        ey[0] = efac*(2.0*q1*q2);
        ey[1] = efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0));
        ey[2] = efac*(2.0*q2);

        // Calculate radial unit vector
        std::vector<T> er = thames::util::vector::mult3<T>(1.0/r, R);

        // Calculate trig of the true longitude
        T cl = thames::util::vector::dot3<T>(er, ex);
        T sl = thames::util::vector::dot3<T>(er, ey);

        // Calculate equinoctial reference frame velocity components
        T hwh = q1*cl - q2*sl;

        // Calculate angular momentum
        std::vector<T> H = thames::util::vector::cross3<T>(R, V);
        T h = thames::util::vector::norm3<T>(H);
        std::vector<T> eh = thames::util::vector::mult3<T>(1.0/h, H);

        // Calculate the effective potential energy
        T ueff = pow(h, 2.0)/(2.0*pow(r, 2.0)) + U;

        // Calculate the generalised angular momentum
        T c = sqrt(2.0*pow(r, 2.0)*ueff);

        // Calculate the generalised semi-latus rectum
        T p = pow(c, 2.0)/mu;

        // Calculate perturbation components
        T Fr = thames::util::vector::dot3<T>(F, er);
        T Fh = thames::util::vector::dot3<T>(F, eh);

        // Calculate non-dimensional quantities
        T zeta = r/p;
        T zetatilde = 1 + zeta;

        // Calculate time derivatives of the second and third elements
        T p1dot = p2*((h - c)/pow(r, 2.0) - r/h*hwh*Fh) + 1.0/c*(r*drdt/c*p1 + zetatilde*p2 + zeta*cl)*(2.0*U - r*Fr) + r/mu*(zeta*p1 + zetatilde*sl)*edot;
        T p2dot = p1*(r/h*hwh*Fh - (h-c)/pow(r, 2.0)) + 1.0/c*(r*drdt/c*p2 - zetatilde*p1 - zeta*sl)*(2.0*U - r*Fr) + r/mu*(zeta*p2 + zetatilde*cl)*edot;

        // Calculate generalised semi-major axis
        T a = pow(mu/pow(nu, 2.0), 1.0/3.0);

        // Calculate alpha
        T alpha = 1.0/(1.0 + sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0)));

        // Calculate time derivative of the generalised mean longitude
        T Ldot = nu + (h - c)/pow(r, 2.0) - r/h*hwh*Fh + (r*drdt*c/pow(mu, 2.0)*zetatilde*alpha)*edot + 1.0/c*(1.0/alpha + alpha*(1.0 - r/a))*(2.0*U - r*Fr);

        // Calculate time derivatives of the remaining elements
        T q1dot = r/(2.0*h)*Fh*(1.0 + pow(q1, 2.0) + pow(q2, 2.0))*sl;
        T q2dot = r/(2.0*h)*Fh*(1.0 + pow(q1, 2.0) + pow(q2, 2.0))*cl;

        // Store derivatives
        geqoedot[0] = nudot;
        geqoedot[1] = p1dot;
        geqoedot[2] = p2dot;
        geqoedot[3] = Ldot;
        geqoedot[4] = q1dot;
        geqoedot[5] = q2dot;
    }
    template void derivative<double>(const std::vector<double>&, std::vector<double>&, const double, const double&, const BasePerturbation<double>&);

    template<class T>
    std::array<T, 6> propagate(T tstart, T tend, T tstep, std::array<T, 6> RV, T mu, const BasePerturbation<T>& perturbation, T atol, T rtol){
        // Transform initial state
        std::array<T, 6> geqoe = thames::conversions::geqoe::cartesian_to_geqoe<T>(tstart, RV, mu, perturbation);

        // Declare derivative function wrapper
        auto derivative_param = [&](const std::array<T, 6>& x, std::array<T, 6>& dxdt, const T time){
            derivative<T>(x, dxdt, time, mu, perturbation);
        };

        // Declare stepper
        boost::numeric::odeint::runge_kutta_cash_karp54<std::array<T, 6>> stepper;
        auto steppercontrolled = boost::numeric::odeint::make_controlled(atol, rtol, stepper);

        // Propagate orbit
        boost::numeric::odeint::integrate_adaptive(steppercontrolled, derivative_param, geqoe, tstart, tend, tstep);

        // Transform final state
        RV = thames::conversions::geqoe::geqoe_to_cartesian<T>(tend, geqoe, mu, perturbation);

        // Return final state
        return RV;
    }
    template std::array<double, 6> propagate<double>(double, double, double, std::array<double, 6>, double, const BasePerturbation<double>&, double, double);

    template<class T>
    std::vector<T> propagate(T tstart, T tend, T tstep, std::vector<T> RV, T mu, const BasePerturbation<T>& perturbation, T atol, T rtol){
        // Transform initial state
        std::vector<T> geqoe = thames::conversions::geqoe::cartesian_to_geqoe<T>(tstart, RV, mu, perturbation);

        // Declare derivative function wrapper
        auto derivative_param = [&](const std::vector<T>& x, std::vector<T>& dxdt, const T time){
            derivative<T>(x, dxdt, time, mu, perturbation);
        };

        // Declare stepper
        boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<T>> stepper;
        auto steppercontrolled = boost::numeric::odeint::make_controlled(atol, rtol, stepper);

        // Propagate orbit
        boost::numeric::odeint::integrate_adaptive(steppercontrolled, derivative_param, geqoe, tstart, tend, tstep);

        // Transform final state
        RV = thames::conversions::geqoe::geqoe_to_cartesian<T>(tend, geqoe, mu, perturbation);

        // Return final state
        return RV;
    }
    template std::vector<double> propagate<double>(double, double, double, std::vector<double>, double, const BasePerturbation<double>&, double, double);

}