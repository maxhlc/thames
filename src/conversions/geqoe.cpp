#include <Eigen/Core>
#include <Eigen/Geometry>

#include "geqoe.h"
#include "keplerian.h"
#include "../perturbations/baseperturbation.h"
#include "../util/root.h"
#include "../types.h"

using namespace thames::types;
using namespace thames::perturbations::baseperturbation;

namespace thames::conversions::geqoe{
    
    Vector6 cartesian_to_geqoe(const double &t, const Vector6 &RV, const double &mu, BasePerturbation &perturbation){
        // Extract position and velocity vectors
        Vector3 R, V;
        R = RV(Eigen::seq(0,2));
        V = RV(Eigen::seq(3,5));

        // Calculate range and range rate
        double r = R.norm();
        double drdt = R.dot(V)/r;

        // Calculate the angular momentum
        Vector3 H = R.cross(V);
        double h = H.norm();

        // Calculate the effective potential energy
        double ueff = pow(h, 2.0)/(2.0*pow(r, 2.0)) + perturbation.potential(t, R);

        // Calculate the total energy
        double e = 0.5*pow(drdt, 2.0) - mu/r + ueff;

        // Calculate the generalised mean motion
        double nu = pow(-2.0*e, 1.5)/mu;

        // Calculate Keplerian elements, and extract angles
        Vector6 keplerian = thames::conversions::keplerian::cartesian_to_keplerian(RV, mu);
        double inc = keplerian[2];
        double raan = keplerian[3];

        // Calculate plane orientation parameters
        double q1 = tan(inc/2.0)*sin(raan);
        double q2 = tan(inc/2.0)*cos(raan);

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

        // Calculate the generalised angular momentum
        double c = sqrt(2.0*pow(r, 2.0)*ueff);

        // Calculate the generalised semi-latus rectum
        double p = pow(c, 2.0)/mu;

        // Calculate remaining non-osculating ellipse parameters
        double pfac1 = (p/r - 1.0);
        double pfac2 = c*drdt/mu;
        double p1 = pfac1*sl - pfac2*cl;
        double p2 = pfac1*cl + pfac2*sl;

        // Calculate generalised semi-major axis and velocity
        double a = pow(mu/pow(nu, 2.0), 1.0/3.0);
        double w = sqrt(mu/a);

        // Calculate generalised mean longitude
        double SCfac1 = mu + c*w - r*pow(drdt, 2.0);
        double SCfac2 = drdt*(c + w*r);
        double S = SCfac1*sl - SCfac2*cl;
        double C = SCfac1*cl + SCfac2*sl;
        double L = atan2(S, C) + (C*p1 - S*p2)/(mu + c*w);

        // Construct GEqOE state vector
        Vector6 geqoe;
        geqoe << nu, p1, p2, L, q1, q2;

        // Return GEqOE state vector
        return geqoe;
    }

    Vector6 geqoe_to_cartesian(const double &t, const Vector6 &geqoe, const double &mu, BasePerturbation &perturbation){
        // Extract elements
        double nu = geqoe[0];
        double p1 = geqoe[1];
        double p2 = geqoe[2];
        double L = geqoe[3];
        double q1 = geqoe[4];
        double q2 = geqoe[5];

        // Calculate generalised eccentric longitude
        auto fk = [p1, p2, L](double k) {return (k + p1*cos(k) - p2*sin(k) - L);};
        auto dfk = [p1, p2, L](double k) {return (1 - p1*sin(k) - p2*cos(k));};
        double k = thames::util::root::newton_raphson(fk, dfk, L);
        double sink = sin(k);
        double cosk = cos(k);

        // Calculate generalised semi-major axis
        double a = pow(mu/pow(nu, 2.0), 1.0/3.0);

        // Calculate range and range rate
        double r = a*(1.0 - p1*sink - p2*cosk);
        double drdt = sqrt(mu*a)/r*(p2*sink - p1*cosk);

        // Calculate trig of the true longitude
        double alpha = 1.0/(1.0 + sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0)));
        double sinl = a/r*(alpha*p1*p2*cosk + (1.0 - alpha*pow(p2, 2.0))*sink - p1);
        double cosl = a/r*(alpha*p1*p2*sink + (1.0 - alpha*pow(p1, 2.0))*cosk - p2);

        // Calculate equinocital reference frame unit vectors
        double efac = 1.0/(1.0 + pow(q1, 2.0) + pow(q2, 2.0));
        Vector3 ex, ey;
        ex << efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0)), efac*(2.0*q1*q2), efac*(-2.0*q1);
        ey << efac*(2.0*q1*q2), efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0)), efac*(2.0*q2);

        // Calculate orbital basis vectors
        Vector3 er = ex*cosl + ey*sinl;
        Vector3 ef = ey*cosl - ex*sinl;

        // Calculate position
        Vector3 R = r*er;

        // Calculate generalised angular momentum
        double c = pow(pow(mu, 2.0)/nu, 1.0/3.0)*sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0));

        // Calculate angular momentum
        double h = sqrt(pow(c, 2.0) - 2.0*pow(r, 2.0)*perturbation.potential(t, R));

        // Calculate velocity
        Vector3 V = drdt*er + h/r*ef;

        // Construct Cartesian state vector
        Vector6 RV;
        RV << R, V;

        // Return Cartesian state vector
        return RV;
    }
}