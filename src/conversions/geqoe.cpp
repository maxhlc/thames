#include "geqoe.h"
#include "keplerian.h"
#include "../perturbations/baseperturbation.h"
#include "../util/root.h"
#include "../util/vector.h"
#include "../types.h"

using namespace thames::types;
using namespace thames::perturbations::baseperturbation;

namespace thames::conversions::geqoe{
    
    template<class real, class vector3, class vector6>
    vector6 cartesian_to_geqoe(const real &t, const vector6 &RV, const real &mu, BasePerturbation<real, vector3> &perturbation){
        // Extract position and velocity vectors
        vector3 R, V;
        R = thames::util::vector::slice<Vector3, Vector6, unsigned int>(RV, 0, 2);
        V = thames::util::vector::slice<Vector3, Vector6, unsigned int>(RV, 3, 5);

        // Calculate range and range rate
        real r = thames::util::vector::norm3<real, vector3>(R);
        real drdt = thames::util::vector::dot3<real, vector3>(R, V)/r;

        // Calculate the angular momentum
        vector3 H = thames::util::vector::cross3<vector3>(R, V);
        real h = thames::util::vector::norm3<real, vector3>(H);

        // Calculate the effective potential energy
        real ueff = pow(h, 2.0)/(2.0*pow(r, 2.0)) + perturbation.potential(t, R);

        // Calculate the total energy
        real e = 0.5*pow(drdt, 2.0) - mu/r + ueff;

        // Calculate the generalised mean motion
        real nu = pow(-2.0*e, 1.5)/mu;

        // Calculate Keplerian elements, and extract angles
        vector6 keplerian = thames::conversions::keplerian::cartesian_to_keplerian<real, vector3, vector6>(RV, mu);
        real inc = keplerian[2];
        real raan = keplerian[3];

        // Calculate plane orientation parameters
        real q1 = tan(inc/2.0)*sin(raan);
        real q2 = tan(inc/2.0)*cos(raan);

        // Calculate equinocital reference frame unit vectors
        real efac = 1.0/(1.0 + pow(q1, 2.0) + pow(q2, 2.0));
        vector3 ex, ey;
        ex[0] = efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0));
        ex[1] = efac*(2.0*q1*q2);
        ex[2] = efac*(-2.0*q1);
        ey[0] = efac*(2.0*q1*q2);
        ey[1] = efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0));
        ey[2] = efac*(2.0*q2);

        // Calculate radial unit vector
        vector3 er = thames::util::vector::mult3<real, vector3>(1.0/r, R);

        // Calculate trig of the true longitude
        real cl = thames::util::vector::dot3<real, vector3>(er, ex);
        real sl = thames::util::vector::dot3<real, vector3>(er, ey);

        // Calculate the generalised angular momentum
        real c = sqrt(2.0*pow(r, 2.0)*ueff);

        // Calculate the generalised semi-latus rectum
        real p = pow(c, 2.0)/mu;

        // Calculate remaining non-osculating ellipse parameters
        real pfac1 = (p/r - 1.0);
        real pfac2 = c*drdt/mu;
        real p1 = pfac1*sl - pfac2*cl;
        real p2 = pfac1*cl + pfac2*sl;

        // Calculate generalised semi-major axis and velocity
        real a = pow(mu/pow(nu, 2.0), 1.0/3.0);
        real w = sqrt(mu/a);

        // Calculate generalised mean longitude
        real SCfac1 = mu + c*w - r*pow(drdt, 2.0);
        real SCfac2 = drdt*(c + w*r);
        real S = SCfac1*sl - SCfac2*cl;
        real C = SCfac1*cl + SCfac2*sl;
        real L = atan2(S, C) + (C*p1 - S*p2)/(mu + c*w);

        // Construct GEqOE state vector
        vector6 geqoe;
        geqoe[0] = nu;
        geqoe[1] = p1;
        geqoe[2] = p2;
        geqoe[3] = L;
        geqoe[4] = q1;
        geqoe[5] = q2;

        // Return GEqOE state vector
        return geqoe;
    }
    template Vector6 cartesian_to_geqoe<double, Vector3, Vector6>(const double&, const Vector6&, const double&, BasePerturbation<double, Vector3>&);

    template<class real, class vector3, class vector6>
    vector6 geqoe_to_cartesian(const real &t, const vector6 &geqoe, const real &mu, BasePerturbation<real, vector3> &perturbation){
        // Extract elements
        real nu = geqoe[0];
        real p1 = geqoe[1];
        real p2 = geqoe[2];
        real L = geqoe[3];
        real q1 = geqoe[4];
        real q2 = geqoe[5];

        // Calculate generalised eccentric longitude
        auto fk = [p1, p2, L](real k) {return (k + p1*cos(k) - p2*sin(k) - L);};
        auto dfk = [p1, p2, L](real k) {return (1 - p1*sin(k) - p2*cos(k));};
        real k = thames::util::root::newton_raphson<real>(fk, dfk, L);
        real sink = sin(k);
        real cosk = cos(k);

        // Calculate generalised semi-major axis
        real a = pow(mu/pow(nu, 2.0), 1.0/3.0);

        // Calculate range and range rate
        real r = a*(1.0 - p1*sink - p2*cosk);
        real drdt = sqrt(mu*a)/r*(p2*sink - p1*cosk);

        // Calculate trig of the true longitude
        real alpha = 1.0/(1.0 + sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0)));
        real sinl = a/r*(alpha*p1*p2*cosk + (1.0 - alpha*pow(p2, 2.0))*sink - p1);
        real cosl = a/r*(alpha*p1*p2*sink + (1.0 - alpha*pow(p1, 2.0))*cosk - p2);

        // Calculate equinocital reference frame unit vectors
        real efac = 1.0/(1.0 + pow(q1, 2.0) + pow(q2, 2.0));
        vector3 ex, ey;
        ex[0] = efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0));
        ex[1] = efac*(2.0*q1*q2);
        ex[2] = efac*(-2.0*q1);
        ey[0] = efac*(2.0*q1*q2);
        ey[1] = efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0));
        ey[2] = efac*(2.0*q2);

        // Calculate orbital basis vectors
        vector3 er, ef;
        for(unsigned int ii=0; ii<3; ii++){
            er[ii] = ex[ii]*cosl + ey[ii]*sinl;
            ef[ii] = ey[ii]*cosl - ex[ii]*sinl;
        }

        // Calculate position
        vector3 R = thames::util::vector::mult3<real, vector3>(r, er);

        // Calculate generalised angular momentum
        real c = pow(pow(mu, 2.0)/nu, 1.0/3.0)*sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0));

        // Calculate angular momentum
        real h = sqrt(pow(c, 2.0) - 2.0*pow(r, 2.0)*perturbation.potential(t, R));

        // Calculate velocity
        vector3 V;
        for(unsigned int ii=0; ii<3; ii++)
            V[ii] = drdt*er[ii] + h/r*ef[ii];

        // Construct Cartesian state vector
        vector6 RV;
        RV[0] = R[0];
        RV[1] = R[1];
        RV[2] = R[2];
        RV[3] = V[0];
        RV[4] = V[1];
        RV[5] = V[2];

        // Return Cartesian state vector
        return RV;
    }
    template Vector6 geqoe_to_cartesian<double, Vector3, Vector6>(const double&, const Vector6&, const double&, BasePerturbation<double, Vector3>&);
    
}