#include <functional>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "state.h"
#include "util.h"
#include "../types.h"
#include "../util/root.h"

using namespace thames::types;

namespace thames::conversions::state{

    Vector6 cartesian_to_keplerian(const Vector6 &RV, const double &mu){
        // TODO: Documentation

        // Extract position and velocity vectors
        Vector3 R, V;
        R = RV(Eigen::seq(0,2));
        V = RV(Eigen::seq(3,5));
        
        // Declare units vectors
        const Vector3 I = {1.0, 0.0, 0.0};
        const Vector3 J = {0.0, 1.0, 0.0};
        const Vector3 K = {0.0, 0.0, 1.0};

        // Set constants
        const double atol = 1e-12;

        // Declare Keplerian elements vector
        Vector6 keplerian;

        // Calculate state magnitudes
        double r = R.norm();
        double v = V.norm();

        // Calculate semi-major axis
        double sma = 1.0/(2.0/r - pow(v, 2.0)/mu);

        // Calculate angular momentum vector and magnitude
        Vector3 H = R.cross(V);
        double h = H.norm();

        // Calculate eccentricity vector and magnitude
        Vector3 E = V.cross(H)/mu - R.normalized();
        double e = E.norm();

        // Calculate inclination
        double inc = acos(K.dot(H)/h);

        // Check for circular and equatorial orbits
        bool e_near = abs(e) < atol;
        bool inc_near = abs(inc) < atol;

        // Calculate right ascension of the ascending node
        Vector3 N = K.cross(H);
        float n = N.norm();
        float raan;
        if (inc_near) {
            raan = 0.0;
        } else {
            raan = acos(I.dot(N)/n);
            if(J.dot(N) < 0.0){
                raan = 2.0*M_PI - raan;
            }
        }

        // Calculate argument of periapsis
        float aop;
        if (inc_near & e_near) {
            aop = 0.0;
        } else if (inc_near) {
            aop = atan2(E[1], E[0]);
            if(H[2] < 0.0){
                aop = 2.0*M_PI - aop;
            }
        } else {
            aop = acos(N.dot(E)/(n*e));
            if(K.dot(E) < 0.0){
                aop = 2.0*M_PI - aop;
            }
        }

        // Calculate true anomaly
        float ta;
        if (inc_near & e_near) {
            ta = acos(R[0]/r);
            if (V[0] > 0.0) {
                ta = 2.0*M_PI - ta;
            }
        } else if (e_near) {
            ta = acos(N.dot(R)/(n*r));
            if (R[2] < 0.0) {
                ta = 2.0*M_PI - ta;
            }
        } else {
            ta = acos(E.dot(R)/(e*r));
            if (R.dot(V) < 0.0) {
                ta = 2.0*M_PI - ta;
            }
        }

        // Store elements in the Keplerian elements vector
        keplerian << sma, e, inc, raan, aop, ta;

        // Return Keplerian elements vector
        return keplerian;
    }

    Vector6 keplerian_to_cartesian(const Vector6 &keplerian, const double &mu){
        // TODO: documentation

        // Extract Keplerian elements
        double sma = keplerian[0];
        double e = keplerian[1];
        double inc = keplerian[2];
        double raan = keplerian[3];
        double aop = keplerian[4];
        double ta = keplerian[5];

        // Calculate true anomaly trigs
        double cta = cos(ta), sta = sin(ta);

        // Calculate eccentric anomaly
        double E = 2.0*atan(sqrt((1.0 - e)/(1.0 + e))*tan(ta/2.0));

        // Calculate radial distance
        double r = sma*(1.0 - pow(e, 2.0))/(1.0 + e*cta);

        // Calculate position and velocity vectors in the orbital frame
        Vector3 o;
        o << r*cta, r*sta, 0.0;

        Vector3 dodt;
        double fac = sqrt(mu*sma)/r;
        dodt << -fac*sin(E), fac*sqrt(1.0 - pow(e, 2.0))*cos(E), 0.0;

        // Calculate rotation matrix
        Matrix33 rot = thames::conversions::util::rot_z(-raan)*
                       thames::conversions::util::rot_x(-inc)*
                       thames::conversions::util::rot_z(-aop);

        // Transform the position and velocity vectors to the inertial frame
        Vector3 R = rot*o;
        Vector3 V = rot*dodt;
        Vector6 RV;
        RV << R, V;

        // Return Cartesian state
        return RV;
    }

    Vector6 cartesian_to_geqoe(const double &t, const Vector6 &RV, const double &mu, const Potential &U){
        // TODO: documentation

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
        double ueff = pow(h, 2.0)/(2.0*pow(r, 2.0)) + U(t, R);

        // Calculate the total energy
        double e = 0.5*pow(drdt, 2.0) - mu/r + ueff;

        // Calculate the generalised mean motion
        double nu = pow(-2.0*e, 1.5)/mu;

        // Calculate Keplerian elements, and extract angles
        Vector6 keplerian = cartesian_to_keplerian(RV, mu);
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

    Vector6 geqoe_to_cartesian(const double &t, const Vector6 &geqoe, const double &mu, const Potential &U){
        // TODO: documentation

        // Extract elements
        double nu = geqoe[0];
        double p1 = geqoe[1];
        double p2 = geqoe[2];
        double L = geqoe[3];
        double q1 = geqoe[4];
        double q2 = geqoe[5];

        // Calculate generalised eccentric longitude
        std::function<double (double)> fk = [p1, p2, L](double k) {return (k + p1*cos(k) - p2*sin(k) - L);};
        double k = thames::util::root::golden_section_search(fk, -M_PI, M_PI);
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
        double h = sqrt(pow(c, 2.0) - 2.0*pow(r, 2.0)*U(t, R));

        // Calculate velocity
        Vector3 V = drdt*er + h/r*ef;

        // Construct Cartesian state vector
        Vector6 RV;
        RV << R, V;

        // Return Cartesian state vector
        return RV;
    }

}