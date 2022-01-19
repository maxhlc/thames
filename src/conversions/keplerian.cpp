#include <array>
#include <cmath>
#include <vector>

#include "keplerian.h"
#include "../util/vector.h"

namespace thames::conversions::keplerian{

    template<class T>
    std::array<T, 6> cartesian_to_keplerian(const std::array<T, 6>& RV, const T& mu){
        // Extract position and velocity vectors
        std::array<T, 3> R, V;
        R = thames::util::vector::slice<T, 6, 3>(RV, 0, 2);
        V = thames::util::vector::slice<T, 6, 3>(RV, 3, 5);
        
        // Declare units vectors
        const std::array<T, 3> I = {1.0, 0.0, 0.0};
        const std::array<T, 3> J = {0.0, 1.0, 0.0};
        const std::array<T, 3> K = {0.0, 0.0, 1.0};

        // Set constants
        const double atol = 1e-12;

        // Declare Keplerian elements vector
        std::array<T, 6> keplerian;

        // Calculate state magnitudes
        T r = thames::util::vector::norm3<T>(R);
        T v = thames::util::vector::norm3<T>(V);

        // Calculate semi-major axis
        T sma = 1.0/(2.0/r - pow(v, 2.0)/mu);

        // Calculate angular momentum vector and magnitude
        std::array<T, 3> H = thames::util::vector::cross3<T>(R, V);
        T h = thames::util::vector::norm3<T>(H);

        // Calculate eccentricity vector and magnitude
        std::array<T, 3> E;
        std::array<T, 3> Etmp = thames::util::vector::cross3<T>(V, H);
        for(unsigned int ii=0; ii<3; ii++)
            E[ii] = Etmp[ii]/mu - R[ii]/r;
        T e = thames::util::vector::norm3<T>(E);

        // Calculate inclination
        T inc = acos(thames::util::vector::dot3<T>(K, H)/h);

        // Check for circular and equatorial orbits
        bool e_near = fabs(e) < atol;
        bool inc_near = fabs(inc) < atol;

        // Calculate right ascension of the ascending node
        std::array<T, 3> N = thames::util::vector::cross3<T>(K, H);
        T n = thames::util::vector::norm3<T>(N);
        T raan;
        if (inc_near) {
            raan = 0.0;
        } else {
            raan = acos(thames::util::vector::dot3<T>(I, N)/n);
            if(thames::util::vector::dot3<T>(J, N) < 0.0){
                raan = 2.0*M_PI - raan;
            }
        }

        // Calculate argument of periapsis
        T aop;
        if (inc_near & e_near) {
            aop = 0.0;
        } else if (inc_near) {
            aop = atan2(E[1], E[0]);
            if(H[2] < 0.0){
                aop = 2.0*M_PI - aop;
            }
        } else {
            aop = acos(thames::util::vector::dot3<T>(N, E)/(n*e));
            if(thames::util::vector::dot3<T>(K, E) < 0.0){
                aop = 2.0*M_PI - aop;
            }
        }

        // Calculate true anomaly
        T ta;
        if (inc_near & e_near) {
            ta = acos(R[0]/r);
            if (V[0] > 0.0) {
                ta = 2.0*M_PI - ta;
            }
        } else if (e_near) {
            ta = acos(thames::util::vector::dot3<T>(N, R)/(n*r));
            if (R[2] < 0.0) {
                ta = 2.0*M_PI - ta;
            }
        } else {
            ta = acos(thames::util::vector::dot3<T>(E, R)/(e*r));
            if (thames::util::vector::dot3<T>(R, V) < 0.0) {
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
    template std::array<double, 6> cartesian_to_keplerian<double>(const std::array<double, 6>&, const double&);

    template<class T>
    std::vector<T> cartesian_to_keplerian(const std::vector<T>& RV, const T& mu){
        // Extract position and velocity vectors
        std::vector<T> R(3), V(3);
        R = thames::util::vector::slice<T>(RV, 0, 2);
        V = thames::util::vector::slice<T>(RV, 3, 5);
        
        // Declare units vectors
        const std::vector<T> I = {1.0, 0.0, 0.0};
        const std::vector<T> J = {0.0, 1.0, 0.0};
        const std::vector<T> K = {0.0, 0.0, 1.0};

        // Set constants
        const double atol = 1e-12;

        // Declare Keplerian elements vector
        std::vector<T> keplerian(6);

        // Calculate state magnitudes
        T r = thames::util::vector::norm3<T>(R);
        T v = thames::util::vector::norm3<T>(V);

        // Calculate semi-major axis
        T sma = 1.0/(2.0/r - pow(v, 2.0)/mu);

        // Calculate angular momentum vector and magnitude
        std::vector<T> H = thames::util::vector::cross3<T>(R, V);
        T h = thames::util::vector::norm3<T>(H);

        // Calculate eccentricity vector and magnitude
        std::vector<T> E(3);
        std::vector<T> Etmp = thames::util::vector::cross3<T>(V, H);
        for(unsigned int ii=0; ii<3; ii++)
            E[ii] = Etmp[ii]/mu - R[ii]/r;
        T e = thames::util::vector::norm3<T>(E);

        // Calculate inclination
        T inc = acos(thames::util::vector::dot3<T>(K, H)/h);

        // Check for circular and equatorial orbits
        bool e_near = fabs(e) < atol;
        bool inc_near = fabs(inc) < atol;

        // Calculate right ascension of the ascending node
        std::vector<T> N = thames::util::vector::cross3<T>(K, H);
        T n = thames::util::vector::norm3<T>(N);
        T raan;
        if (inc_near) {
            raan = 0.0;
        } else {
            raan = acos(thames::util::vector::dot3<T>(I, N)/n);
            if(thames::util::vector::dot3<T>(J, N) < 0.0){
                raan = 2.0*M_PI - raan;
            }
        }

        // Calculate argument of periapsis
        T aop;
        if (inc_near & e_near) {
            aop = 0.0;
        } else if (inc_near) {
            aop = atan2(E[1], E[0]);
            if(H[2] < 0.0){
                aop = 2.0*M_PI - aop;
            }
        } else {
            aop = acos(thames::util::vector::dot3<T>(N, E)/(n*e));
            if(thames::util::vector::dot3<T>(K, E) < 0.0){
                aop = 2.0*M_PI - aop;
            }
        }

        // Calculate true anomaly
        T ta;
        if (inc_near & e_near) {
            ta = acos(R[0]/r);
            if (V[0] > 0.0) {
                ta = 2.0*M_PI - ta;
            }
        } else if (e_near) {
            ta = acos(thames::util::vector::dot3<T>(N, R)/(n*r));
            if (R[2] < 0.0) {
                ta = 2.0*M_PI - ta;
            }
        } else {
            ta = acos(thames::util::vector::dot3<T>(E, R)/(e*r));
            if (thames::util::vector::dot3<T>(R, V) < 0.0) {
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
    template std::vector<double> cartesian_to_keplerian<double>(const std::vector<double>&, const double&);

    template<class T>
    std::array<T, 6> keplerian_to_cartesian(const std::array<T, 6>& keplerian, const T& mu){
        // Extract Keplerian elements
        T sma = keplerian[0];
        T e = keplerian[1];
        T inc = keplerian[2];
        T raan = keplerian[3];
        T aop = keplerian[4];
        T ta = keplerian[5];

        // Calculate angle trigs
        T cinc = cos(inc), sinc = sin(inc);
        T craan = cos(raan), sraan = sin(raan);
        T caop = cos(aop), saop = sin(aop);
        T cta = cos(ta), sta = sin(ta);

        // Calculate eccentric anomaly
        T E = 2.0*atan(sqrt((1.0 - e)/(1.0 + e))*tan(ta/2.0));

        // Calculate radial distance
        T r = sma*(1.0 - pow(e, 2.0))/(1.0 + e*cta);

        // Calculate position and velocity vectors in the orbital frame
        std::array<T, 3> o;
        o[0] = r*cta;
        o[1] = r*sta;
        o[2] = 0.0;

        std::array<T, 3> dodt;
        T fac = sqrt(mu*sma)/r;
        dodt[0] = -fac*sin(E);
        dodt[1] = fac*sqrt(1.0 - pow(e, 2.0))*cos(E);
        dodt[2] = 0.0;

        // Calculate rotation angles
        T ang[3][2];
        ang[0][0] = caop*craan - saop*cinc*sraan;
        ang[0][1] = -saop*craan - caop*cinc*sraan;
        ang[1][0] = caop*sraan + saop*cinc*craan;
        ang[1][1] = caop*cinc*craan - saop*sraan;
        ang[2][0] = saop*sinc;
        ang[2][1] = caop*sinc;

        // Transform the position and velocity vectors to the inertial frame
        std::array<T, 6> RV;
        for(unsigned int ii=0; ii<3; ii++){
            RV[ii] = o[0]*ang[ii][0] + o[1]*ang[ii][1];
            RV[ii+3] = dodt[0]*ang[ii][0] + dodt[1]*ang[ii][1];
        }

        // Return Cartesian state
        return RV;
    }
    template std::array<double, 6> keplerian_to_cartesian<double>(const std::array<double, 6>&, const double&);

    template<class T>
    std::vector<T> keplerian_to_cartesian(const std::vector<T>& keplerian, const T& mu){
        // Extract Keplerian elements
        T sma = keplerian[0];
        T e = keplerian[1];
        T inc = keplerian[2];
        T raan = keplerian[3];
        T aop = keplerian[4];
        T ta = keplerian[5];

        // Calculate angle trigs
        T cinc = cos(inc), sinc = sin(inc);
        T craan = cos(raan), sraan = sin(raan);
        T caop = cos(aop), saop = sin(aop);
        T cta = cos(ta), sta = sin(ta);

        // Calculate eccentric anomaly
        T E = 2.0*atan(sqrt((1.0 - e)/(1.0 + e))*tan(ta/2.0));

        // Calculate radial distance
        T r = sma*(1.0 - pow(e, 2.0))/(1.0 + e*cta);

        // Calculate position and velocity vectors in the orbital frame
        std::vector<T> o(3);
        o[0] = r*cta;
        o[1] = r*sta;
        o[2] = 0.0;

        std::vector<T> dodt(3);
        T fac = sqrt(mu*sma)/r;
        dodt[0] = -fac*sin(E);
        dodt[1] = fac*sqrt(1.0 - pow(e, 2.0))*cos(E);
        dodt[2] = 0.0;

        // Calculate rotation angles
        T ang[3][2];
        ang[0][0] = caop*craan - saop*cinc*sraan;
        ang[0][1] = -saop*craan - caop*cinc*sraan;
        ang[1][0] = caop*sraan + saop*cinc*craan;
        ang[1][1] = caop*cinc*craan - saop*sraan;
        ang[2][0] = saop*sinc;
        ang[2][1] = caop*sinc;

        // Transform the position and velocity vectors to the inertial frame
        std::vector<T> RV(6);
        for(unsigned int ii=0; ii<3; ii++){
            RV[ii] = o[0]*ang[ii][0] + o[1]*ang[ii][1];
            RV[ii+3] = dodt[0]*ang[ii][0] + dodt[1]*ang[ii][1];
        }

        // Return Cartesian state
        return RV;
    }
    template std::vector<double> keplerian_to_cartesian<double>(const std::vector<double>&, const double&);

}