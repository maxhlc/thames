#include <cmath>

#include "keplerian.h"
#include "../util/vector.h"
#include "../types.h"

using namespace thames::types;

namespace thames::conversions::keplerian{

    template<class real, class vector3, class vector6>
    vector6 cartesian_to_keplerian(const vector6 &RV, const real &mu){
        // Extract position and velocity vectors
        vector3 R, V;
        R = thames::util::vector::slice<Vector3, Vector6, unsigned int>(RV, 0, 2);
        V = thames::util::vector::slice<Vector3, Vector6, unsigned int>(RV, 3, 5);
        
        // Declare units vectors
        const vector3 I = {1.0, 0.0, 0.0};
        const vector3 J = {0.0, 1.0, 0.0};
        const vector3 K = {0.0, 0.0, 1.0};

        // Set constants
        const real atol = 1e-12;

        // Declare Keplerian elements vector
        vector6 keplerian;

        // Calculate state magnitudes
        real r = thames::util::vector::norm3<real, vector3>(R);
        real v = thames::util::vector::norm3<real, vector3>(V);

        // Calculate semi-major axis
        real sma = 1.0/(2.0/r - pow(v, 2.0)/mu);

        // Calculate angular momentum vector and magnitude
        vector3 H = thames::util::vector::cross3<vector3>(R, V);
        real h = thames::util::vector::norm3<real, vector3>(H);

        // Calculate eccentricity vector and magnitude
        vector3 E;
        vector3 Etmp = thames::util::vector::cross3<vector3>(V, H);
        for(unsigned int ii=0; ii<3; ii++)
            E[ii] = Etmp[ii]/mu - R[ii]/r;
        real e = thames::util::vector::norm3<real, vector3>(E);

        // Calculate inclination
        real inc = acos(thames::util::vector::dot3<real, vector3>(K, H)/h);

        // Check for circular and equatorial orbits
        bool e_near = fabs(e) < atol;
        bool inc_near = fabs(inc) < atol;

        // Calculate right ascension of the ascending node
        vector3 N = thames::util::vector::cross3<vector3>(K, H);
        real n = thames::util::vector::norm3<real, vector3>(N);
        real raan;
        if (inc_near) {
            raan = 0.0;
        } else {
            raan = acos(thames::util::vector::dot3<real, vector3>(I, N)/n);
            if(thames::util::vector::dot3<real, vector3>(J, N) < 0.0){
                raan = 2.0*M_PI - raan;
            }
        }

        // Calculate argument of periapsis
        real aop;
        if (inc_near & e_near) {
            aop = 0.0;
        } else if (inc_near) {
            aop = atan2(E[1], E[0]);
            if(H[2] < 0.0){
                aop = 2.0*M_PI - aop;
            }
        } else {
            aop = acos(thames::util::vector::dot3<real, vector3>(N, E)/(n*e));
            if(thames::util::vector::dot3<real, vector3>(K, E) < 0.0){
                aop = 2.0*M_PI - aop;
            }
        }

        // Calculate true anomaly
        real ta;
        if (inc_near & e_near) {
            ta = acos(R[0]/r);
            if (V[0] > 0.0) {
                ta = 2.0*M_PI - ta;
            }
        } else if (e_near) {
            ta = acos(thames::util::vector::dot3<real, vector3>(N, R)/(n*r));
            if (R[2] < 0.0) {
                ta = 2.0*M_PI - ta;
            }
        } else {
            ta = acos(thames::util::vector::dot3<real, vector3>(E, R)/(e*r));
            if (thames::util::vector::dot3<real, vector3>(R, V) < 0.0) {
                ta = 2.0*M_PI - ta;
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
    template Vector6 cartesian_to_keplerian<double, Vector3, Vector6>(const Vector6&, const double&);

    template<class real, class vector3, class vector6>
    vector6 keplerian_to_cartesian(const vector6 &keplerian, const real &mu){
        // Extract Keplerian elements
        real sma = keplerian[0];
        real e = keplerian[1];
        real inc = keplerian[2];
        real raan = keplerian[3];
        real aop = keplerian[4];
        real ta = keplerian[5];

        // Calculate angle trigs
        real cinc = cos(inc), sinc = sin(inc);
        real craan = cos(raan), sraan = sin(raan);
        real caop = cos(aop), saop = sin(aop);
        real cta = cos(ta), sta = sin(ta);

        // Calculate eccentric anomaly
        real E = 2.0*atan(sqrt((1.0 - e)/(1.0 + e))*tan(ta/2.0));

        // Calculate radial distance
        real r = sma*(1.0 - pow(e, 2.0))/(1.0 + e*cta);

        // Calculate position and velocity vectors in the orbital frame
        vector3 o;
        o[0] = r*cta;
        o[1] = r*sta;
        o[2] = 0.0;

        vector3 dodt;
        real fac = sqrt(mu*sma)/r;
        dodt[0] = -fac*sin(E);
        dodt[1] = fac*sqrt(1.0 - pow(e, 2.0))*cos(E);
        dodt[2] = 0.0;

        // Calculate rotation angles
        real ang[3][2];
        ang[0][0] = caop*craan - saop*cinc*sraan;
        ang[0][1] = -saop*craan - caop*cinc*sraan;
        ang[1][0] = caop*sraan + saop*cinc*craan;
        ang[1][1] = caop*cinc*craan - saop*sraan;
        ang[2][0] = saop*sinc;
        ang[2][1] = caop*sinc;

        // Transform the position and velocity vectors to the inertial frame
        vector6 RV;
        for(unsigned int ii=0; ii<3; ii++){
            RV[ii] = o[0]*ang[ii][0] + o[1]*ang[ii][1];
            RV[ii+3] = dodt[0]*ang[ii][0] + dodt[1]*ang[ii][1];
        }

        // Return Cartesian state
        return RV;
    }
    template Vector6 keplerian_to_cartesian<double, Vector3, Vector6>(const Vector6&, const double&);

}