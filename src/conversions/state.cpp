#include <array>
#include <iostream>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>

#include "state.h"
#include "../types.h"

using namespace thames::types;

namespace thames::conversions::state{

    Vector6 cartesian_to_keplerian(const Vector3 &R, const Vector3 &V, const double &mu){
        // TODO: Documentation
        
        // Declare units vectors
        const Vector3 I = {1.0f, 0.0f, 0.0f};
        const Vector3 J = {0.0f, 1.0f, 0.0f};
        const Vector3 K = {0.0f, 0.0f, 1.0f};

        // Set constants
        const double atol = 1e-12;

        // Declare Keplerian elements vector
        Vector6 keplerian;

        // Calculate state magnitudes
        double r = R.norm();
        double v = V.norm();

        // Calculate semi-major axis
        double sma = 1.0f/(2.0f/r - pow(v, 2.0f)/mu);

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
            raan = 0.0f;
        } else {
            raan = acos(I.dot(N)/n);
            if(J.dot(N) < 0.0f){
                raan = 2*M_PI - raan;
            }
        }

        // Calculate argument of periapsis
        float aop;
        if (inc_near & e_near) {
            aop = 0.0f;
        } else if (inc_near) {
            aop = atan2(E[1], E[0]);
            if(H[2] < 0.0f){
                aop = 2*M_PI - aop;
            }
        } else {
            aop = acos(N.dot(E)/(n*e));
            if(K.dot(E) < 0.0){
                aop = 2*M_PI - aop;
            }
        }

        // Calculate true anomaly
        float ta;
        if (inc_near & e_near) {
            ta = acos(R[0]/r);
            if (V[0] > 0.0f) {
                ta = 2*M_PI - ta;
            }
        } else if (e_near) {
            ta = acos(N.dot(R)/(n*r));
            if (R[2] < 0.0f) {
                ta = 2*M_PI - ta;
            }
        } else {
            ta = acos(E.dot(R)/(e*r));
            if (R.dot(V) < 0.0f) {
                ta = 2*M_PI - ta;
            }
        }

        // Store elements in the Keplerian elements vector
        keplerian[0] = sma;
        keplerian[1] = e;
        keplerian[2] = inc;
        keplerian[3] = raan;
        keplerian[4] = aop;
        keplerian[5] = ta;

        // Return Keplerian elements vector
        return keplerian;
    }

}