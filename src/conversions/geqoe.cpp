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
#include <functional>
#include <memory>
#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../../include/conversions/geqoe.h"
#include "../../include/conversions/keplerian.h"
#include "../../include/perturbations/baseperturbation.h"
#include "../../include/util/root.h"
#include "../../include/vector/arithmeticoverloads.h"
#include "../../include/vector/geometry.h"

namespace thames::conversions::geqoe{

    using namespace thames::perturbations::baseperturbation;
    using namespace thames::vector::arithmeticoverloads;

    ////////////
    // Arrays //
    ////////////

    template<class T>
    std::array<T, 6> cartesian_to_geqoe(const T& t, const std::array<T, 6>& RV, const T& mu, const std::shared_ptr<const BasePerturbation<T>> perturbation){
        // Extract position and velocity vectors
        std::array<T, 3> R = {RV[0], RV[1], RV[2]};
        std::array<T, 3> V = {RV[3], RV[4], RV[5]};

        // Calculate range and range rate
        T r = thames::vector::geometry::norm3(R);
        T drdt = thames::vector::geometry::dot3(R, V)/r;

        // Calculate the angular momentum
        std::array<T, 3> H = thames::vector::geometry::cross3(R, V);
        T h = thames::vector::geometry::norm3(H);

        // Calculate the effective potential energy
        T ueff = pow(h, 2.0)/(2.0*pow(r, 2.0)) + perturbation->potential(t, R);

        // Calculate the total energy
        T e = 0.5*pow(drdt, 2.0) - mu/r + ueff;

        // Calculate the generalised mean motion
        T nu = pow(-2.0*e, 1.5)/mu;

        // // Calculate plane orientation parameters
        T q1 = H[0]/(h + H[2]);
        T q2 = -H[1]/(h + H[2]);

        // Calculate equinocital reference frame unit vectors
        T efac = 1.0/(1.0 + pow(q1, 2.0) + pow(q2, 2.0));
        std::array<T, 3> ex = {
            efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0)),
            efac*(2.0*q1*q2),
            efac*(-2.0*q1)
        };
        std::array<T, 3> ey = {
            efac*(2.0*q1*q2),
            efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0)),
            efac*(2.0*q2)
        };

        // Calculate radial unit vector
        std::array<T, 3> er = R/r;

        // Calculate trig of the true longitude
        T cl = thames::vector::geometry::dot3(er, ex);
        T sl = thames::vector::geometry::dot3(er, ey);

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
        std::array<T, 6> geqoe = {
            nu,
            p1,
            p2,
            L,
            q1,
            q2
        };

        // Return GEqOE state vector
        return geqoe;
    }
    template std::array<double, 6> cartesian_to_geqoe<double>(const double&, const std::array<double, 6>&, const double&, const std::shared_ptr<const BasePerturbation<double>> perturbation);

    template<class T>
    std::array<T, 6> geqoe_to_cartesian(const T& t, const std::array<T, 6>& geqoe, const T& mu, const std::shared_ptr<const BasePerturbation<T>> perturbation){
        // Extract elements
        T nu = geqoe[0];
        T p1 = geqoe[1];
        T p2 = geqoe[2];
        T L = geqoe[3];
        T q1 = geqoe[4];
        T q2 = geqoe[5];

        // Calculate generalised eccentric longitude
        std::function<double (double)> fk = [p1, p2, L](T k) {return (k + p1*cos(k) - p2*sin(k) - L);};
        std::function<double (double)> dfk = [p1, p2, L](T k) {return (1 - p1*sin(k) - p2*cos(k));};
        T k = thames::util::root::newton_raphson(fk, dfk, L);
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
        std::array<T, 3> ex = {
            efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0)),
            efac*(2.0*q1*q2),
            efac*(-2.0*q1)
        };
        std::array<T, 3> ey = {
            efac*(2.0*q1*q2),
            efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0)),
            efac*(2.0*q2)
        };

        // Calculate orbital basis vectors
        std::array<T, 3> er = ex*cosl + ey*sinl;
        std::array<T, 3> ef = ey*cosl - ex*sinl;

        // Calculate position
        std::array<T, 3> R = r*er;

        // Calculate generalised angular momentum
        T c = pow(pow(mu, 2.0)/nu, 1.0/3.0)*sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0));

        // Calculate angular momentum
        T h = sqrt(pow(c, 2.0) - 2.0*pow(r, 2.0)*perturbation->potential(t, R));

        // Calculate velocity
        std::array<T, 3> V = drdt*er + h/r*ef;

        // Construct Cartesian state vector
        std::array<T, 6> RV = {
            R[0],
            R[1],
            R[2],
            V[0],
            V[1],
            V[2]
        };

        // Return Cartesian state vector
        return RV;
    }
    template std::array<double, 6> geqoe_to_cartesian<double>(const double&, const std::array<double, 6>&, const double&, const std::shared_ptr<const BasePerturbation<double>> perturbation);

    /////////////
    // Vectors //
    /////////////

    template<class T>
    std::vector<T> cartesian_to_geqoe(const T& t, const std::vector<T>& RV, const T& mu, const std::shared_ptr<const BasePerturbation<T>> perturbation){
        // Extract position and velocity vectors
        std::vector<T> R = {RV[0], RV[1], RV[2]};
        std::vector<T> V = {RV[3], RV[4], RV[5]};

        // Calculate range and range rate
        T r = thames::vector::geometry::norm3(R);
        T drdt = thames::vector::geometry::dot3(R, V)/r;

        // Calculate the angular momentum
        std::vector<T> H = thames::vector::geometry::cross3(R, V);
        T h = thames::vector::geometry::norm3(H);

        // Calculate the effective potential energy
        T ueff = pow(h, 2.0)/(2.0*pow(r, 2.0)) + perturbation->potential(t, R);

        // Calculate the total energy
        T e = 0.5*pow(drdt, 2.0) - mu/r + ueff;

        // Calculate the generalised mean motion
        T nu = pow(-2.0*e, 1.5)/mu;

        // Calculate plane orientation parameters
        T q1 = H[0]/(h + H[2]);
        T q2 = -H[1]/(h + H[2]);

        // Calculate equinocital reference frame unit vectors
        T efac = 1.0/(1.0 + pow(q1, 2.0) + pow(q2, 2.0));
        std::vector<T> ex = {
            efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0)),
            efac*(2.0*q1*q2),
            efac*(-2.0*q1)
        };
        std::vector<T> ey = {
            efac*(2.0*q1*q2),
            efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0)),
            efac*(2.0*q2)
        };

        // Calculate radial unit vector
        std::vector<T> er = R/r;

        // Calculate trig of the true longitude
        T cl = thames::vector::geometry::dot3(er, ex);
        T sl = thames::vector::geometry::dot3(er, ey);

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
        std::vector<T> geqoe = {
            nu,
            p1,
            p2,
            L,
            q1,
            q2
        };

        // Return GEqOE state vector
        return geqoe;
    }
    template std::vector<double> cartesian_to_geqoe<double>(const double& t, const std::vector<double>& RV, const double& mu, const std::shared_ptr<const BasePerturbation<double>> perturbation);

    template<class T>
    std::vector<T> geqoe_to_cartesian(const T& t, const std::vector<T>& geqoe, const T& mu, const std::shared_ptr<const BasePerturbation<T>> perturbation){
        // Extract elements
        T nu = geqoe[0];
        T p1 = geqoe[1];
        T p2 = geqoe[2];
        T L = geqoe[3];
        T q1 = geqoe[4];
        T q2 = geqoe[5];

        // Calculate generalised eccentric longitude
        std::function<double (double)> fk = [p1, p2, L](T k) {return (k + p1*cos(k) - p2*sin(k) - L);};
        std::function<double (double)> dfk = [p1, p2, L](T k) {return (1 - p1*sin(k) - p2*cos(k));};
        T k = thames::util::root::newton_raphson(fk, dfk, L);
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
        std::vector<T> ex = {
            efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0)),
            efac*(2.0*q1*q2),
            efac*(-2.0*q1)
        };
        std::vector<T> ey = {
            efac*(2.0*q1*q2),
            efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0)),
            efac*(2.0*q2)
        };

        // Calculate orbital basis vectors
        std::vector<T> er = ex*cosl + ey*sinl;
        std::vector<T> ef = ey*cosl - ex*sinl;

        // Calculate position
        std::vector<T> R = r*er;

        // Calculate generalised angular momentum
        T c = pow(pow(mu, 2.0)/nu, 1.0/3.0)*sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0));

        // Calculate angular momentum
        T h = sqrt(pow(c, 2.0) - 2.0*pow(r, 2.0)*perturbation->potential(t, R));

        // Calculate velocity
        std::vector<T> V = drdt*er + h/r*ef;

        // Construct Cartesian state vector
        std::vector<T> RV = {
            R[0],
            R[1],
            R[2],
            V[0],
            V[1],
            V[2]
        };

        // Return Cartesian state vector
        return RV;
    }
    template std::vector<double> geqoe_to_cartesian<double>(const double&, const std::vector<double>&, const double&, const std::shared_ptr<const BasePerturbation<double>> perturbation);

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::polynomial;

    template<class T, template<class> class P>
    std::vector<P<T>> cartesian_to_geqoe(const T& t, const std::vector<P<T>>& RV, const T& mu, const std::shared_ptr<const BasePerturbationPolynomial<T, P>> perturbation){
        // Extract position and velocity vectors
        std::vector<P<T>> R = {RV[0], RV[1], RV[2]};
        std::vector<P<T>> V = {RV[3], RV[4], RV[5]};

        // Calculate range and range rate
        P<T> r = thames::vector::geometry::norm3(R);
        P<T> drdt = thames::vector::geometry::dot3(R, V)/r;

        // Calculate the angular momentum
        std::vector<P<T>> H = thames::vector::geometry::cross3(R, V);
        P<T> h = thames::vector::geometry::norm3(H);

        // Calculate the effective potential energy
        P<T> ueff = pow(h, 2.0)/(2.0*pow(r, 2.0)) + perturbation->potential(t, R);

        // Calculate the total energy
        P<T> e = 0.5*pow(drdt, 2.0) - mu/r + ueff;

        // Calculate the generalised mean motion
        P<T> nu = pow(-2.0*e, 1.5)/mu;

        // Calculate plane orientation parameters
        P<T> q1 = H[0]/(h + H[2]);
        P<T> q2 = -H[1]/(h + H[2]);

        // Calculate equinocital reference frame unit vectors
        P<T> efac = 1.0/(1.0 + pow(q1, 2.0) + pow(q2, 2.0));
        std::vector<P<T>> ex = {
            efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0)),
            efac*(2.0*q1*q2),
            efac*(-2.0*q1)
        };
        std::vector<P<T>> ey = {
            efac*(2.0*q1*q2),
            efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0)),
            efac*(2.0*q2)
        };

        // Calculate radial unit vector
        std::vector<P<T>> er = R/r;

        // Calculate trig of the true longitude
        P<T> cl = thames::vector::geometry::dot3(er, ex);
        P<T> sl = thames::vector::geometry::dot3(er, ey);

        // Calculate the generalised angular momentum
        P<T> c = sqrt(2.0*pow(r, 2.0)*ueff);

        // Calculate the generalised semi-latus rectum
        P<T> p = pow(c, 2.0)/mu;

        // Calculate remaining non-osculating ellipse parameters
        P<T> pfac1 = (p/r - 1.0);
        P<T> pfac2 = c*drdt/mu;
        P<T> p1 = pfac1*sl - pfac2*cl;
        P<T> p2 = pfac1*cl + pfac2*sl;

        // Calculate generalised semi-major axis and velocity
        P<T> a = pow(mu/pow(nu, 2.0), 1.0/3.0);
        P<T> w = sqrt(mu/a);

        // Calculate generalised mean longitude
        P<T> SCfac1 = mu + c*w - r*pow(drdt, 2.0);
        P<T> SCfac2 = drdt*(c + w*r);
        P<T> S = SCfac1*sl - SCfac2*cl;
        P<T> C = SCfac1*cl + SCfac2*sl;
        P<T> L = atan2(S, C) + (C*p1 - S*p2)/(mu + c*w);

        // Construct GEqOE state vector
        std::vector<P<T>> geqoe = {
            nu,
            p1,
            p2,
            L,
            q1,
            q2
        };

        // Return GEqOE state vector
        return geqoe;
    }
    template std::vector<taylor_polynomial<double>> cartesian_to_geqoe(const double& t, const std::vector<taylor_polynomial<double>>& RV, const double& mu, const std::shared_ptr<const BasePerturbationPolynomial<double, taylor_polynomial>> perturbation);
    template std::vector<chebyshev_polynomial<double>> cartesian_to_geqoe(const double& t, const std::vector<chebyshev_polynomial<double>>& RV, const double& mu, const std::shared_ptr<const BasePerturbationPolynomial<double, chebyshev_polynomial>> perturbation);

    template<class T, template<class> class P>
    std::vector<P<T>> geqoe_to_cartesian(const T& t, const std::vector<P<T>>& geqoe, const T& mu, const std::shared_ptr<const BasePerturbationPolynomial<T, P>> perturbation){
        // Extract elements
        P<T> nu = geqoe[0];
        P<T> p1 = geqoe[1];
        P<T> p2 = geqoe[2];
        P<T> L = geqoe[3];
        P<T> q1 = geqoe[4];
        P<T> q2 = geqoe[5];

        // Calculate generalised eccentric longitude
        std::function<P<T> (P<T>)> fk = [p1, p2, L](P<T> k) {return (k + p1*cos(k) - p2*sin(k) - L);};
        std::function<P<T> (P<T>)> dfk = [p1, p2, L](P<T> k) {return (1 - p1*sin(k) - p2*cos(k));};
        P<T> k = thames::util::root::newton_raphson(fk, dfk, L);
        P<T> sink = sin(k);
        P<T> cosk = cos(k);

        // Calculate generalised semi-major axis
        P<T> a = pow(mu/pow(nu, 2.0), 1.0/3.0);

        // Calculate range and range rate
        P<T> r = a*(1.0 - p1*sink - p2*cosk);
        P<T> drdt = sqrt(mu*a)/r*(p2*sink - p1*cosk);

        // Calculate trig of the true longitude
        P<T> alpha = 1.0/(1.0 + sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0)));
        P<T> sinl = a/r*(alpha*p1*p2*cosk + (1.0 - alpha*pow(p2, 2.0))*sink - p1);
        P<T> cosl = a/r*(alpha*p1*p2*sink + (1.0 - alpha*pow(p1, 2.0))*cosk - p2);

        // Calculate equinocital reference frame unit vectors
        P<T> efac = 1.0/(1.0 + pow(q1, 2.0) + pow(q2, 2.0));
        std::vector<P<T>> ex = {
            efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0)),
            efac*(2.0*q1*q2),
            efac*(-2.0*q1)
        };
        std::vector<P<T>> ey = {
            efac*(2.0*q1*q2),
            efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0)),
            efac*(2.0*q2)
        };

        // Calculate orbital basis vectors
        std::vector<P<T>> er = ex*cosl + ey*sinl;
        std::vector<P<T>> ef = ey*cosl - ex*sinl;

        // Calculate position
        std::vector<P<T>> R = r*er;

        // Calculate generalised angular momentum
        P<T> c = pow(pow(mu, 2.0)/nu, 1.0/3.0)*sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0));

        // Calculate angular momentum
        P<T> h = sqrt(pow(c, 2.0) - 2.0*pow(r, 2.0)*perturbation->potential(t, R));

        // Calculate velocity
        std::vector<P<T>> V = drdt*er + h/r*ef;

        // Construct Cartesian state vector
        std::vector<P<T>> RV = {
            R[0],
            R[1],
            R[2],
            V[0],
            V[1],
            V[2]
        };

        // Return Cartesian state vector
        return RV;
    }
    template std::vector<taylor_polynomial<double>> geqoe_to_cartesian(const double&, const std::vector<taylor_polynomial<double>>&, const double&, const std::shared_ptr<const BasePerturbationPolynomial<double, taylor_polynomial>> perturbation);
    template std::vector<chebyshev_polynomial<double>> geqoe_to_cartesian(const double&, const std::vector<chebyshev_polynomial<double>>&, const double&, const std::shared_ptr<const BasePerturbationPolynomial<double, chebyshev_polynomial>> perturbation);

    #endif

}