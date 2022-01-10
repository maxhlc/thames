#include <Eigen/Core>
#include <Eigen/Geometry>

#include "keplerian.h"
#include "util.h"
#include "../types.h"

using namespace thames::types;

namespace thames::conversions::keplerian{

    Vector6 cartesian_to_keplerian(const Vector6 &RV, const double &mu){
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
        double n = N.norm();
        double raan;
        if (inc_near) {
            raan = 0.0;
        } else {
            raan = acos(I.dot(N)/n);
            if(J.dot(N) < 0.0){
                raan = 2.0*M_PI - raan;
            }
        }

        // Calculate argument of periapsis
        double aop;
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
        double ta;
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

}