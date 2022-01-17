#include <array>
#include <cmath>

#include "geqoe.h"
#include "keplerian.h"
#include "../perturbations/baseperturbation.h"
#include "../util/root.h"
#include "../util/vector.h"

using namespace thames::perturbations::baseperturbation;

namespace thames::conversions::geqoe{
    
    template<class T>
    std::array<T, 6> cartesian_to_geqoe(const T& t, const std::array<T, 6>& RV, const T& mu, const BasePerturbation<T>& perturbation){
        // Extract position and velocity vectors
        std::array<T, 3> R, V;
        R = thames::util::vector::slice<T, 6, 3>(RV, 0, 2);
        V = thames::util::vector::slice<T, 6, 3>(RV, 3, 5);

        // Calculate range and range rate
        T r = thames::util::vector::norm3<T>(R);
        T drdt = thames::util::vector::dot3<T>(R, V)/r;

        // Calculate the angular momentum
        std::array<T, 3> H = thames::util::vector::cross3<T>(R, V);
        T h = thames::util::vector::norm3<T>(H);

        // Calculate the effective potential energy
        T ueff = pow(h, 2.0)/(2.0*pow(r, 2.0)) + perturbation.potential(t, R);

        // Calculate the total energy
        T e = 0.5*pow(drdt, 2.0) - mu/r + ueff;

        // Calculate the generalised mean motion
        T nu = pow(-2.0*e, 1.5)/mu;

        // Calculate Keplerian elements, and extract angles
        std::array<T, 6> keplerian = thames::conversions::keplerian::cartesian_to_keplerian<T>(RV, mu);
        T inc = keplerian[2];
        T raan = keplerian[3];

        // Calculate plane orientation parameters
        T q1 = tan(inc/2.0)*sin(raan);
        T q2 = tan(inc/2.0)*cos(raan);

        // Calculate equinocital reference frame unit vectors
        T efac = 1.0/(1.0 + pow(q1, 2.0) + pow(q2, 2.0));
        std::array<T, 3> ex, ey;
        ex[0] = efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0));
        ex[1] = efac*(2.0*q1*q2);
        ex[2] = efac*(-2.0*q1);
        ey[0] = efac*(2.0*q1*q2);
        ey[1] = efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0));
        ey[2] = efac*(2.0*q2);

        // Calculate radial unit vector
        std::array<T, 3> er = thames::util::vector::mult3<T>(1.0/r, R);

        // Calculate trig of the true longitude
        T cl = thames::util::vector::dot3<T>(er, ex);
        T sl = thames::util::vector::dot3<T>(er, ey);

        // Calculate the generalised angular momentum
        T c = sqrt(2.0*pow(r, 2.0)*ueff);

        // Calculate the generalised semi-latus rectum
        T p = pow(c, 2.0)/mu;

        // Calculate remaining non-osculating ellipse parameters
        T pfac1 = (p/r - 1.0);
        T pfac2 = c*drdt/mu;
        T p1 = pfac1*sl - pfac2*cl;
        T p2 = pfac1*cl + pfac2*sl;

        // Calculate generalised semi-major axis and velocity
        T a = pow(mu/pow(nu, 2.0), 1.0/3.0);
        T w = sqrt(mu/a);

        // Calculate generalised mean longitude
        T SCfac1 = mu + c*w - r*pow(drdt, 2.0);
        T SCfac2 = drdt*(c + w*r);
        T S = SCfac1*sl - SCfac2*cl;
        T C = SCfac1*cl + SCfac2*sl;
        T L = atan2(S, C) + (C*p1 - S*p2)/(mu + c*w);

        // Construct GEqOE state vector
        std::array<T, 6> geqoe;
        geqoe[0] = nu;
        geqoe[1] = p1;
        geqoe[2] = p2;
        geqoe[3] = L;
        geqoe[4] = q1;
        geqoe[5] = q2;

        // Return GEqOE state vector
        return geqoe;
    }
    template std::array<double, 6> cartesian_to_geqoe<double>(const double&, const std::array<double, 6>&, const double&, const BasePerturbation<double>&);

    template<class T>
    std::array<T, 6> geqoe_to_cartesian(const T& t, const std::array<T, 6>& geqoe, const T& mu, const BasePerturbation<T>& perturbation){
        // Extract elements
        T nu = geqoe[0];
        T p1 = geqoe[1];
        T p2 = geqoe[2];
        T L = geqoe[3];
        T q1 = geqoe[4];
        T q2 = geqoe[5];

        // Calculate generalised eccentric longitude
        auto fk = [p1, p2, L](T k) {return (k + p1*cos(k) - p2*sin(k) - L);};
        auto dfk = [p1, p2, L](T k) {return (1 - p1*sin(k) - p2*cos(k));};
        T k = thames::util::root::newton_raphson<T>(fk, dfk, L);
        T sink = sin(k);
        T cosk = cos(k);

        // Calculate generalised semi-major axis
        T a = pow(mu/pow(nu, 2.0), 1.0/3.0);

        // Calculate range and range rate
        T r = a*(1.0 - p1*sink - p2*cosk);
        T drdt = sqrt(mu*a)/r*(p2*sink - p1*cosk);

        // Calculate trig of the true longitude
        T alpha = 1.0/(1.0 + sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0)));
        T sinl = a/r*(alpha*p1*p2*cosk + (1.0 - alpha*pow(p2, 2.0))*sink - p1);
        T cosl = a/r*(alpha*p1*p2*sink + (1.0 - alpha*pow(p1, 2.0))*cosk - p2);

        // Calculate equinocital reference frame unit vectors
        T efac = 1.0/(1.0 + pow(q1, 2.0) + pow(q2, 2.0));
        std::array<T, 3> ex, ey;
        ex[0] = efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0));
        ex[1] = efac*(2.0*q1*q2);
        ex[2] = efac*(-2.0*q1);
        ey[0] = efac*(2.0*q1*q2);
        ey[1] = efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0));
        ey[2] = efac*(2.0*q2);

        // Calculate orbital basis vectors
        std::array<T, 3> er, ef;
        for(unsigned int ii=0; ii<3; ii++){
            er[ii] = ex[ii]*cosl + ey[ii]*sinl;
            ef[ii] = ey[ii]*cosl - ex[ii]*sinl;
        }

        // Calculate position
        std::array<T, 3> R = thames::util::vector::mult3<T>(r, er);

        // Calculate generalised angular momentum
        T c = pow(pow(mu, 2.0)/nu, 1.0/3.0)*sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0));

        // Calculate angular momentum
        T h = sqrt(pow(c, 2.0) - 2.0*pow(r, 2.0)*perturbation.potential(t, R));

        // Calculate velocity
        std::array<T, 3> V;
        for(unsigned int ii=0; ii<3; ii++)
            V[ii] = drdt*er[ii] + h/r*ef[ii];

        // Construct Cartesian state vector
        std::array<T, 6> RV;
        RV[0] = R[0];
        RV[1] = R[1];
        RV[2] = R[2];
        RV[3] = V[0];
        RV[4] = V[1];
        RV[5] = V[2];

        // Return Cartesian state vector
        return RV;
    }
    template std::array<double, 6> geqoe_to_cartesian<double>(const double&, const std::array<double, 6>&, const double&, const BasePerturbation<double>&);
    
}