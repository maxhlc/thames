#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

#include "geqoe.h"
#include "../types.h"
#include "../conversions/geqoe.h"
#include "../perturbations/baseperturbation.h"

using namespace thames::types;
using namespace thames::perturbations::baseperturbation;

namespace thames::propagators::geqoe{

    void derivative(const Vector6 &geqoe, Vector6 &geqoedot, const double t, const double &mu, BasePerturbation &perturbation){
        // Extract elements
        double nu = geqoe[0];
        double p1 = geqoe[1];
        double p2 = geqoe[2];
        double q1 = geqoe[4];
        double q2 = geqoe[5];

        // Calculate Cartesian state and extract vectors
        Vector6 RV = thames::conversions::geqoe::geqoe_to_cartesian(t, geqoe, mu, perturbation);
        Vector3 R, V;
        R = RV(Eigen::seq(0,2));
        V = RV(Eigen::seq(3,5));

        // Calculate range and range rate
        double r = R.norm();
        double drdt = R.dot(V)/r;

        // Calculate perturbations
        double U = perturbation.potential(t, R);
        double Ut = perturbation.potential_derivative(t, R, V);
        Vector3 F = perturbation.acceleration_total(t, R, V);
        Vector3 P = perturbation.acceleration_nonpotential(t, R, V);

        // Calculate time derivative of total energy
        double edot = Ut + P.dot(V);

        // Calculate time derivative of nu
        double nudot = -3.0*pow(nu/pow(mu, 2.0), 1.0/3.0)*edot;

        // Calculate equinocital reference frame unit vectors
        double efac = 1.0/(1.0 + pow(q1, 2.0) + pow(q2, 2.0));
        Vector3 ex, ey;
        ex << efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0)), efac*(2.0*q1*q2), efac*(-2.0*q1);
        ey << efac*(2.0*q1*q2), efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0)), efac*(2.0*q2);

        // Calculate radial unit vector
        Vector3 er = R.normalized();

        // Calculate trig of the true longitude
        double cl = er.dot(ex);
        double sl = er.dot(ey);

        // Calculate equinoctial reference frame velocity components
        double hwh = q1*cl - q2*sl;

        // Calculate angular momentum
        Vector3 H = R.cross(V);
        double h = H.norm();
        Vector3 eh = H.normalized();

        // Calculate the effective potential energy
        double ueff = pow(h, 2.0)/(2.0*pow(r, 2.0)) + U;

        // Calculate the generalised angular momentum
        double c = sqrt(2.0*pow(r, 2.0)*ueff);

        // Calculate the generalised semi-latus rectum
        double p = pow(c, 2.0)/mu;

        // Calculate perturbation components
        double Fr = F.dot(er);
        double Fh = F.dot(eh);

        // Calculate non-dimensional quantities
        double zeta = r/p;
        double zetatilde = 1 + zeta;

        // Calculate time derivatives of the second and third elements
        double p1dot = p2*((h - c)/pow(r, 2.0) - r/h*hwh*Fh) + 1.0/c*(r*drdt/c*p1 + zetatilde*p2 + zeta*cl)*(2.0*U - r*Fr) + r/mu*(zeta*p1 + zetatilde*sl)*edot;
        double p2dot = p1*(r/h*hwh*Fh - (h-c)/pow(r, 2.0)) + 1.0/c*(r*drdt/c*p2 - zetatilde*p1 - zeta*sl)*(2.0*U - r*Fr) + r/mu*(zeta*p2 + zetatilde*cl)*edot;

        // Calculate generalised semi-major axis
        double a = pow(mu/pow(nu, 2.0), 1.0/3.0);

        // Calculate alpha
        double alpha = 1.0/(1.0 + sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0)));

        // Calculate time derivative of the generalised mean longitude
        double Ldot = nu + (h - c)/pow(r, 2.0) - r/h*hwh*Fh + (r*drdt*c/pow(mu, 2.0)*zetatilde*alpha)*edot + 1.0/c*(1.0/alpha + alpha*(1.0 - r/a))*(2.0*U - r*Fr);

        // Calculate time derivatives of the remaining elements
        double q1dot = r/(2.0*h)*Fh*(1.0 + pow(q1, 2.0) + pow(q2, 2.0))*sl;
        double q2dot = r/(2.0*h)*Fh*(1.0 + pow(q1, 2.0) + pow(q2, 2.0))*cl;

        // Store derivatives
        geqoedot << nudot, p1dot, p2dot, Ldot, q1dot, q2dot;
    }

    Vector6 propagate(double tstart, double tend, double tstep, Vector6 RV, double mu, BasePerturbation &perturbation, double atol, double rtol){
        // Transform initial state
        Vector6 geqoe = thames::conversions::geqoe::cartesian_to_geqoe(tstart, RV, mu, perturbation);

        // Declare derivative function wrapper
        auto derivative_param = [&](const Vector6 &x, Vector6 &dxdt, const double time){
            derivative(x, dxdt, time, mu, perturbation);
        };

        // Declare stepper
        boost::numeric::odeint::runge_kutta_cash_karp54<Vector6> stepper;
        auto steppercontrolled = boost::numeric::odeint::make_controlled(atol, rtol, stepper);

        // Propagate orbit
        boost::numeric::odeint::integrate_adaptive(steppercontrolled, derivative_param, geqoe, tstart, tend, tstep);

        // Transform final state
        RV = thames::conversions::geqoe::geqoe_to_cartesian(tend, geqoe, mu, perturbation);

        // Return final state
        return RV;
    }

}