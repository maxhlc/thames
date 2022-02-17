/*
MIT License

Copyright (c) 2021-2022 Max Hallgarten La Casta

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <array>
#include <cmath>
#include <vector>

#include "keplerian.h"
#include "../vector/arithmeticoverloads.h"
#include "../vector/geometry.h"

using namespace thames::vector::arithmeticoverloads;

namespace thames::conversions::keplerian{

    ////////////
    // Arrays //
    ////////////

    template<class T>
    std::array<T, 6> cartesian_to_keplerian(const std::array<T, 6>& RV, const T& mu){
        // Extract position and velocity vectors
        std::array<T, 3> R = {RV[0], RV[1], RV[2]};
        std::array<T, 3> V = {RV[3], RV[4], RV[5]};
        
        // Declare units vectors
        const std::array<T, 3> K = {0.0, 0.0, 1.0};

        // Set constants
        const double atol = 1e-12;

        // Calculate state magnitudes
        T r = thames::vector::geometry::norm3(R);
        T v = thames::vector::geometry::norm3(V);

        // Calculate semi-major axis
        T sma = 1.0/(2.0/r - pow(v, 2.0)/mu);

        // Calculate angular momentum vector and magnitude
        std::array<T, 3> H = thames::vector::geometry::cross3(R, V);
        T h = thames::vector::geometry::norm3(H);

        // Calculate eccentricity vector and magnitude
        std::array<T, 3> E = thames::vector::geometry::cross3(V, H)/mu - R/r;
        T e = thames::vector::geometry::norm3(E);

        // Calculate inclination
        T inc = acos(H[2]/h);

        // Check for circular and equatorial orbits
        bool e_near = fabs(e) < atol;
        bool inc_near = fabs(inc) < atol;

        // Calculate right ascension of the ascending node
        std::array<T, 3> N = thames::vector::geometry::cross3(K, H);
        T n = thames::vector::geometry::norm3(N);
        T raan;
        if (inc_near) {
            raan = 0.0;
        } else {
            raan = acos(N[0]/n);
            if(N[1] < 0.0){
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
            aop = acos(thames::vector::geometry::dot3(N, E)/(n*e));
            if(E[2] < 0.0){
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
            ta = acos(thames::vector::geometry::dot3(N, R)/(n*r));
            if (R[2] < 0.0) {
                ta = 2.0*M_PI - ta;
            }
        } else {
            ta = acos(thames::vector::geometry::dot3(E, R)/(e*r));
            if (thames::vector::geometry::dot3(R, V) < 0.0) {
                ta = 2.0*M_PI - ta;
            }
        }

        // Store elements in the Keplerian elements vector
        std::array<T, 6> keplerian = {
            sma,
            e,
            inc,
            raan,
            aop,
            ta
        };

        // Return Keplerian elements vector
        return keplerian;
    }
    template std::array<double, 6> cartesian_to_keplerian<double>(const std::array<double, 6>&, const double&);

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
        std::array<T, 3> o = {
            r*cta,
            r*sta,
            0.0
        };

        T fac = sqrt(mu*sma)/r;
        std::array<T, 3> dodt = {
            -fac*sin(E),
            fac*sqrt(1.0 - pow(e, 2.0))*cos(E),
            0.0
        };

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

    /////////////
    // Vectors //
    /////////////

    template<class T>
    std::vector<T> cartesian_to_keplerian(const std::vector<T>& RV, const T& mu){
        // Extract position and velocity vectors
        std::vector<T> R = {RV[0], RV[1], RV[2]};
        std::vector<T> V = {RV[3], RV[4], RV[5]};
        
        // Declare units vectors
        const std::vector<T> K = {0.0, 0.0, 1.0};

        // Set constants
        const double atol = 1e-12;

        // Calculate state magnitudes
        T r = thames::vector::geometry::norm3(R);
        T v = thames::vector::geometry::norm3(V);

        // Calculate semi-major axis
        T sma = 1.0/(2.0/r - pow(v, 2.0)/mu);

        // Calculate angular momentum vector and magnitude
        std::vector<T> H = thames::vector::geometry::cross3(R, V);
        T h = thames::vector::geometry::norm3(H);

        // Calculate eccentricity vector and magnitude
        std::vector<T> E = thames::vector::geometry::cross3(V, H)/mu - R/r;
        T e = thames::vector::geometry::norm3(E);

        // Calculate inclination
        T inc = acos(H[2]/h);

        // Check for circular and equatorial orbits
        bool e_near = fabs(e) < atol;
        bool inc_near = fabs(inc) < atol;

        // Calculate right ascension of the ascending node
        std::vector<T> N = thames::vector::geometry::cross3(K, H);
        T n = thames::vector::geometry::norm3(N);
        T raan;
        if (inc_near) {
            raan = 0.0;
        } else {
            raan = acos(N[0]/n);
            if(N[1] < 0.0){
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
            aop = acos(thames::vector::geometry::dot3(N, E)/(n*e));
            if(E[2] < 0.0){
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
            ta = acos(thames::vector::geometry::dot3(N, R)/(n*r));
            if (R[2] < 0.0) {
                ta = 2.0*M_PI - ta;
            }
        } else {
            ta = acos(thames::vector::geometry::dot3(E, R)/(e*r));
            if (thames::vector::geometry::dot3(R, V) < 0.0) {
                ta = 2.0*M_PI - ta;
            }
        }

        // Store elements in the Keplerian elements vector
        std::vector<T> keplerian = {
            sma,
            e,
            inc,
            raan,
            aop,
            ta
        };

        // Return Keplerian elements vector
        return keplerian;
    }
    template std::vector<double> cartesian_to_keplerian<double>(const std::vector<double>&, const double&);

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
        std::vector<T> o = {
            r*cta,
            r*sta,
            0.0
        };

        T fac = sqrt(mu*sma)/r;
        std::vector<T> dodt = {
            -fac*sin(E),
            fac*sqrt(1.0 - pow(e, 2.0))*cos(E),
            0.0
        };

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