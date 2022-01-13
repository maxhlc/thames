#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

#include "geqoe.h"
#include "../conversions/geqoe.h"
#include "../perturbations/baseperturbation.h"
#include "../util/vector.h"
#include "../types.h"

using namespace thames::types;
using namespace thames::perturbations::baseperturbation;

namespace thames::propagators::geqoe{

    template<class real, class vector3, class vector6>
    void derivative(const vector6 &geqoe, vector6 &geqoedot, const real t, const real &mu, BasePerturbation<real, vector3> &perturbation){
        // Extract elements
        real nu = geqoe[0];
        real p1 = geqoe[1];
        real p2 = geqoe[2];
        real q1 = geqoe[4];
        real q2 = geqoe[5];

        // Calculate Cartesian state and extract vectors
        vector6 RV = thames::conversions::geqoe::geqoe_to_cartesian<real, vector3, vector6>(t, geqoe, mu, perturbation);
        vector3 R, V;
        R = thames::util::vector::slice<vector3, vector6, unsigned int>(RV, 0, 2);
        V = thames::util::vector::slice<vector3, vector6, unsigned int>(RV, 3, 5);

        // Calculate range and range rate
        real r = thames::util::vector::norm3<real, vector3>(R);
        real drdt = thames::util::vector::dot3<real, vector3>(R, V)/r;

        // Calculate perturbations
        real U = perturbation.potential(t, R);
        real Ut = perturbation.potential_derivative(t, R, V);
        vector3 F = perturbation.acceleration_total(t, R, V);
        vector3 P = perturbation.acceleration_nonpotential(t, R, V);

        // Calculate time derivative of total energy
        real edot = Ut + thames::util::vector::dot3<real, vector3>(P, V);

        // Calculate time derivative of nu
        real nudot = -3.0*pow(nu/pow(mu, 2.0), 1.0/3.0)*edot;

        // Calculate equinocital reference frame unit vectors
        real efac = 1.0/(1.0 + pow(q1, 2.0) + pow(q2, 2.0));
        vector3 ex, ey;
        ex[0] = efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0));
        ex[1] = efac*(2.0*q1*q2);
        ex[2] = efac*(-2.0*q1);
        ey[0] = efac*(2.0*q1*q2);
        ey[1] = efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0));
        ey[2] = efac*(2.0*q2);

        // Calculate radial unit vector
        vector3 er = thames::util::vector::mult3<real, vector3>(1.0/r, R);

        // Calculate trig of the true longitude
        real cl = thames::util::vector::dot3<real, vector3>(er, ex);
        real sl = thames::util::vector::dot3<real, vector3>(er, ey);

        // Calculate equinoctial reference frame velocity components
        real hwh = q1*cl - q2*sl;

        // Calculate angular momentum
        vector3 H = thames::util::vector::cross3<vector3>(R, V);
        real h = thames::util::vector::norm3<real, vector3>(H);
        vector3 eh = thames::util::vector::mult3<real, vector3>(1.0/h, H);

        // Calculate the effective potential energy
        real ueff = pow(h, 2.0)/(2.0*pow(r, 2.0)) + U;

        // Calculate the generalised angular momentum
        real c = sqrt(2.0*pow(r, 2.0)*ueff);

        // Calculate the generalised semi-latus rectum
        real p = pow(c, 2.0)/mu;

        // Calculate perturbation components
        real Fr = thames::util::vector::dot3<real, vector3>(F, er);
        real Fh = thames::util::vector::dot3<real, vector3>(F, eh);

        // Calculate non-dimensional quantities
        real zeta = r/p;
        real zetatilde = 1 + zeta;

        // Calculate time derivatives of the second and third elements
        real p1dot = p2*((h - c)/pow(r, 2.0) - r/h*hwh*Fh) + 1.0/c*(r*drdt/c*p1 + zetatilde*p2 + zeta*cl)*(2.0*U - r*Fr) + r/mu*(zeta*p1 + zetatilde*sl)*edot;
        real p2dot = p1*(r/h*hwh*Fh - (h-c)/pow(r, 2.0)) + 1.0/c*(r*drdt/c*p2 - zetatilde*p1 - zeta*sl)*(2.0*U - r*Fr) + r/mu*(zeta*p2 + zetatilde*cl)*edot;

        // Calculate generalised semi-major axis
        real a = pow(mu/pow(nu, 2.0), 1.0/3.0);

        // Calculate alpha
        real alpha = 1.0/(1.0 + sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0)));

        // Calculate time derivative of the generalised mean longitude
        real Ldot = nu + (h - c)/pow(r, 2.0) - r/h*hwh*Fh + (r*drdt*c/pow(mu, 2.0)*zetatilde*alpha)*edot + 1.0/c*(1.0/alpha + alpha*(1.0 - r/a))*(2.0*U - r*Fr);

        // Calculate time derivatives of the remaining elements
        real q1dot = r/(2.0*h)*Fh*(1.0 + pow(q1, 2.0) + pow(q2, 2.0))*sl;
        real q2dot = r/(2.0*h)*Fh*(1.0 + pow(q1, 2.0) + pow(q2, 2.0))*cl;

        // Store derivatives
        geqoedot[0] = nudot;
        geqoedot[1] = p1dot;
        geqoedot[2] = p2dot;
        geqoedot[3] = Ldot;
        geqoedot[4] = q1dot;
        geqoedot[5] = q2dot;
    }
    template void derivative(const Vector6&, Vector6&, const double, const double&, BasePerturbation<double, Vector3>&);

    template<class real, class vector3, class vector6>
    vector6 propagate(real tstart, real tend, real tstep, vector6 RV, real mu, BasePerturbation<real, vector3> &perturbation, real atol, real rtol){
        // Transform initial state
        vector6 geqoe = thames::conversions::geqoe::cartesian_to_geqoe<real, vector3, vector6>(tstart, RV, mu, perturbation);

        // Declare derivative function wrapper
        auto derivative_param = [&](const vector6 &x, vector6 &dxdt, const real time){
            derivative<real, vector3, vector6>(x, dxdt, time, mu, perturbation);
        };

        // Declare stepper
        boost::numeric::odeint::runge_kutta_cash_karp54<vector6> stepper;
        auto steppercontrolled = boost::numeric::odeint::make_controlled(atol, rtol, stepper);

        // Propagate orbit
        boost::numeric::odeint::integrate_adaptive(steppercontrolled, derivative_param, geqoe, tstart, tend, tstep);

        // Transform final state
        RV = thames::conversions::geqoe::geqoe_to_cartesian<real, vector3, vector6>(tend, geqoe, mu, perturbation);

        // Return final state
        return RV;
    }
    template Vector6 propagate<double, Vector3, Vector6>(double, double, double, Vector6, double, BasePerturbation<double, Vector3>&, double, double);

}