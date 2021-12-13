#include "geqoe.h"
#include "../types.h"
#include "../conversions/state.h"

using namespace thames::types;

namespace thames::propagators::geqoe{

    void derivative(const Vector6 &geqoe, Vector6 &geqoedot, const double t, const double &mu, Potential &U_func, PotentialDerivative &Ut_func, Force &F_func, Force &P_func){
        // TODO: documentation
        // TODO: optimisations

        // Extract elements
        double nu = geqoe[0];
        double p1 = geqoe[1];
        double p2 = geqoe[2];
        double L = geqoe[3];
        double q1 = geqoe[4];
        double q2 = geqoe[5];

        // Calculate Cartesian state and extract vectors
        Vector6 RV = thames::conversions::state::geqoe_to_cartesian(t, geqoe, mu, U_func);
        Vector3 R, V;
        R = RV(Eigen::seq(0,2));
        V = RV(Eigen::seq(3,5));

        // Calculate range and range rate
        double r = R.norm();
        double drdt = R.dot(V)/r;

        // Calculate perturbations
        double U = U_func(t, R);
        double Ut = Ut_func(t, R, V);
        Vector3 F = F_func(t, R, V);
        Vector3 P = P_func(t, R, V);

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
        double hwx = cl;
        double hwy = sl;
        double hwh = q1*cl - q2*sl;

        // Calculate angular momentum
        Vector3 H = R.cross(V);
        double h = H.norm();
        Vector3 eh = H.normalized();

        // Calculate the effective potential energy
        double ueff = pow(h, 2.0)/(2.0*pow(r, 2.0)) + U_func(t, R);

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

}