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

#include <boost/numeric/odeint.hpp>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Integrators/rk4.h"
#include "../../external/smart-uq/include/Integrators/rk45.h"
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../../include/propagators/geqoe.h"
#include "../../include/propagators/options.h"
#include "../../include/conversions/geqoe.h"
#include "../../include/perturbations/baseperturbation.h"
#include "../../include/util/root.h"
#include "../../include/vector/arithmeticoverloads.h"
#include "../../include/vector/geometry.h"

namespace thames::propagators {

    using thames::constants::statetypes::GEQOE;
    using thames::perturbations::baseperturbation::BasePerturbation;
    using namespace thames::vector::arithmeticoverloads;

    template<class T>
    GEqOEPropagator<T>::GEqOEPropagator(const T& mu, BasePerturbation<T>* const perturbation, const DimensionalFactors<T>* factors) : BasePropagator<T>(mu, perturbation, factors, GEQOE) {

    }

    ////////////
    // Arrays //
    ////////////

    template<class T>
    void GEqOEPropagator<T>::derivative(const std::array<T, 6>& geqoe, std::array<T, 6>& geqoedot, const T t) const {
        // Calculate factors
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;

        // Extract elements
        T nu = geqoe[0];
        T p1 = geqoe[1];
        T p2 = geqoe[2];
        T L = geqoe[3];
        T q1 = geqoe[4];
        T q2 = geqoe[5];

        // Calculate generalised eccentric longitude
        std::function<T (T)> fk = [p1, p2, L](T k) {return (k + p1*cos(k) - p2*sin(k) - L);};
        std::function<T (T)> dfk = [p1, p2, L](T k) {return (1 - p1*sin(k) - p2*cos(k));};
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
        T h = sqrt(pow(c, 2.0) - 2.0*pow(r, 2.0)*m_perturbation->potential(t, R));

        // Calculate velocity
        std::array<T, 3> V = drdt*er + h/r*ef;

        // Calculate perturbations
        T U = m_perturbation->potential(t, R);
        T Ut = m_perturbation->potential_derivative(t, R, V);
        std::array<T, 3> F = m_perturbation->acceleration_total(t, R, V);
        std::array<T, 3> P = m_perturbation->acceleration_nonpotential(t, R, V);

        // Calculate time derivative of total energy
        T edot = Ut + thames::vector::geometry::dot3(P, V);

        // Calculate time derivative of nu
        T nudot = -3.0*pow(nu/pow(mu, 2.0), 1.0/3.0)*edot;

        // Calculate trig of the true longitude
        T cl = thames::vector::geometry::dot3(er, ex);
        T sl = thames::vector::geometry::dot3(er, ey);

        // Calculate equinoctial reference frame velocity components
        T hwh = q1*cl - q2*sl;

        // Calculate angular momentum
        std::array<T, 3> H = thames::vector::geometry::cross3(R, V);
        std::array<T, 3> eh = H/h;

        // Calculate the generalised semi-latus rectum
        T p = pow(c, 2.0)/mu;

        // Calculate perturbation components
        T Fr = thames::vector::geometry::dot3(F, er);
        T Fh = thames::vector::geometry::dot3(F, eh);

        // Calculate non-dimensional quantities
        T zeta = r/p;
        T zetatilde = 1 + zeta;

        // Calculate time derivatives of the second and third elements
        T p1dot = p2*((h - c)/pow(r, 2.0) - r/h*hwh*Fh) + 1.0/c*(r*drdt/c*p1 + zetatilde*p2 + zeta*cl)*(2.0*U - r*Fr) + r/mu*(zeta*p1 + zetatilde*sl)*edot;
        T p2dot = p1*(r/h*hwh*Fh - (h-c)/pow(r, 2.0)) + 1.0/c*(r*drdt/c*p2 - zetatilde*p1 - zeta*sl)*(2.0*U - r*Fr) + r/mu*(zeta*p2 + zetatilde*cl)*edot;

        // Calculate time derivative of the generalised mean longitude
        T Ldot = nu + (h - c)/pow(r, 2.0) - r/h*hwh*Fh + (r*drdt*c/pow(mu, 2.0)*zetatilde*alpha)*edot + 1.0/c*(1.0/alpha + alpha*(1.0 - r/a))*(2.0*U - r*Fr);

        // Calculate time derivatives of the remaining elements
        T q1dot = r/(2.0*h)*Fh*(1.0 + pow(q1, 2.0) + pow(q2, 2.0))*sl;
        T q2dot = r/(2.0*h)*Fh*(1.0 + pow(q1, 2.0) + pow(q2, 2.0))*cl;

        // Store derivatives
        geqoedot = {
            nudot,
            p1dot,
            p2dot,
            Ldot,
            q1dot,
            q2dot
        };
    }

    /////////////
    // Vectors //
    /////////////

    template<class T>
    void GEqOEPropagator<T>::derivative(const std::vector<T>& geqoe, std::vector<T>& geqoedot, const T t) const {
        // Calculate factors
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;

        // Extract elements
        T nu = geqoe[0];
        T p1 = geqoe[1];
        T p2 = geqoe[2];
        T L = geqoe[3];
        T q1 = geqoe[4];
        T q2 = geqoe[5];

        // Calculate generalised eccentric longitude
        std::function<T (T)> fk = [p1, p2, L](T k) {return (k + p1*cos(k) - p2*sin(k) - L);};
        std::function<T (T)> dfk = [p1, p2, L](T k) {return (1 - p1*sin(k) - p2*cos(k));};
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
        T h = sqrt(pow(c, 2.0) - 2.0*pow(r, 2.0)*m_perturbation->potential(t, R));

        // Calculate velocity
        std::vector<T> V = drdt*er + h/r*ef;

        // Calculate perturbations
        T U = m_perturbation->potential(t, R);
        T Ut = m_perturbation->potential_derivative(t, R, V);
        std::vector<T> F = m_perturbation->acceleration_total(t, R, V);
        std::vector<T> P = m_perturbation->acceleration_nonpotential(t, R, V);

        // Calculate time derivative of total energy
        T edot = Ut + thames::vector::geometry::dot3(P, V);

        // Calculate time derivative of nu
        T nudot = -3.0*pow(nu/pow(mu, 2.0), 1.0/3.0)*edot;

        // Calculate trig of the true longitude
        T cl = thames::vector::geometry::dot3(er, ex);
        T sl = thames::vector::geometry::dot3(er, ey);

        // Calculate equinoctial reference frame velocity components
        T hwh = q1*cl - q2*sl;

        // Calculate angular momentum
        std::vector<T> H = thames::vector::geometry::cross3(R, V);
        std::vector<T> eh = H/h;

        // Calculate the generalised semi-latus rectum
        T p = pow(c, 2.0)/mu;

        // Calculate perturbation components
        T Fr = thames::vector::geometry::dot3(F, er);
        T Fh = thames::vector::geometry::dot3(F, eh);

        // Calculate non-dimensional quantities
        T zeta = r/p;
        T zetatilde = 1 + zeta;

        // Calculate time derivatives of the second and third elements
        T p1dot = p2*((h - c)/pow(r, 2.0) - r/h*hwh*Fh) + 1.0/c*(r*drdt/c*p1 + zetatilde*p2 + zeta*cl)*(2.0*U - r*Fr) + r/mu*(zeta*p1 + zetatilde*sl)*edot;
        T p2dot = p1*(r/h*hwh*Fh - (h-c)/pow(r, 2.0)) + 1.0/c*(r*drdt/c*p2 - zetatilde*p1 - zeta*sl)*(2.0*U - r*Fr) + r/mu*(zeta*p2 + zetatilde*cl)*edot;

        // Calculate time derivative of the generalised mean longitude
        T Ldot = nu + (h - c)/pow(r, 2.0) - r/h*hwh*Fh + (r*drdt*c/pow(mu, 2.0)*zetatilde*alpha)*edot + 1.0/c*(1.0/alpha + alpha*(1.0 - r/a))*(2.0*U - r*Fr);

        // Calculate time derivatives of the remaining elements
        T q1dot = r/(2.0*h)*Fh*(1.0 + pow(q1, 2.0) + pow(q2, 2.0))*sl;
        T q2dot = r/(2.0*h)*Fh*(1.0 + pow(q1, 2.0) + pow(q2, 2.0))*cl;

        // Store derivatives
        geqoedot = {
            nudot,
            p1dot,
            p2dot,
            Ldot,
            q1dot,
            q2dot
        };
    }

    template class GEqOEPropagator<double>;

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::integrator;
    using namespace smartuq::polynomial;
    using thames::perturbations::baseperturbation::BasePerturbationPolynomial;

    template<class T, template<class> class P>
    GEqOEPropagatorPolynomialDynamics<T, P>::GEqOEPropagatorPolynomialDynamics(const T& mu, BasePerturbationPolynomial<T, P>* const perturbation, const DimensionalFactors<T>* factors) : BasePropagatorPolynomialDynamics<T, P>("GEqOE", mu, perturbation, factors) {

    }

    template<class T, template<class> class P>
    GEqOEPropagatorPolynomialDynamics<T, P>::~GEqOEPropagatorPolynomialDynamics() {

    }

    template<class T, template<class> class W>
    int GEqOEPropagatorPolynomialDynamics<T, W>::evaluate(const T& t, const std::vector<W<T>>& geqoe, std::vector<W<T>>& geqoedot) const {
        // Calculate factors
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;

        // Extract elements
        W<T> nu = geqoe[0];
        W<T> p1 = geqoe[1];
        W<T> p2 = geqoe[2];
        W<T> L = geqoe[3];
        W<T> q1 = geqoe[4];
        W<T> q2 = geqoe[5];

        // Calculate generalised eccentric longitude
        std::function<W<T> (W<T>)> fk = [p1, p2, L](W<T> k) {return (k + p1*cos(k) - p2*sin(k) - L);};
        std::function<W<T> (W<T>)> dfk = [p1, p2, L](W<T> k) {return (1 - p1*sin(k) - p2*cos(k));};
        W<T> k = thames::util::root::newton_raphson(fk, dfk, L);
        W<T> sink = sin(k);
        W<T> cosk = cos(k);

        // Calculate generalised semi-major axis
        W<T> a = pow(mu/pow(nu, 2.0), 1.0/3.0);

        // Calculate range and range rate
        W<T> r = a*(1.0 - p1*sink - p2*cosk);
        W<T> drdt = sqrt(mu*a)/r*(p2*sink - p1*cosk);

        // Calculate trig of the true longitude
        W<T> alpha = 1.0/(1.0 + sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0)));
        W<T> sinl = a/r*(alpha*p1*p2*cosk + (1.0 - alpha*pow(p2, 2.0))*sink - p1);
        W<T> cosl = a/r*(alpha*p1*p2*sink + (1.0 - alpha*pow(p1, 2.0))*cosk - p2);

        // Calculate equinocital reference frame unit vectors
        W<T> efac = 1.0/(1.0 + pow(q1, 2.0) + pow(q2, 2.0));
        std::vector<W<T>> ex = {
            efac*(1.0 - pow(q1, 2.0) + pow(q2, 2.0)),
            efac*(2.0*q1*q2),
            efac*(-2.0*q1)
        };
        std::vector<W<T>> ey = {
            efac*(2.0*q1*q2),
            efac*(1.0 + pow(q1, 2.0) - pow(q2, 2.0)),
            efac*(2.0*q2)
        };

        // Calculate orbital basis vectors
        std::vector<W<T>> er = ex*cosl + ey*sinl;
        std::vector<W<T>> ef = ey*cosl - ex*sinl;

        // Calculate position
        std::vector<W<T>> R = r*er;

        // Calculate generalised angular momentum
        W<T> c = pow(pow(mu, 2.0)/nu, 1.0/3.0)*sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0));

        // Calculate angular momentum
        W<T> h = sqrt(pow(c, 2.0) - 2.0*pow(r, 2.0)*m_perturbation->potential(t, R));

        // Calculate velocity
        std::vector<W<T>> V = drdt*er + h/r*ef;

        // Calculate perturbations
        W<T> U = m_perturbation->potential(t, R);
        W<T> Ut = m_perturbation->potential_derivative(t, R, V);
        std::vector<W<T>> F = m_perturbation->acceleration_total(t, R, V);
        std::vector<W<T>> P = m_perturbation->acceleration_nonpotential(t, R, V);

        // Calculate time derivative of total energy
        W<T> edot = Ut + thames::vector::geometry::dot3(P, V);

        // Calculate time derivative of nu
        W<T> nudot = -3.0*pow(nu/pow(mu, 2.0), 1.0/3.0)*edot;

        // Calculate trig of the true longitude
        W<T> cl = thames::vector::geometry::dot3(er, ex);
        W<T> sl = thames::vector::geometry::dot3(er, ey);

        // Calculate equinoctial reference frame velocity components
        W<T> hwh = q1*cl - q2*sl;

        // Calculate angular momentum
        std::vector<W<T>> H = thames::vector::geometry::cross3(R, V);
        std::vector<W<T>> eh = H/h;

        // Calculate the generalised semi-latus rectum
        W<T> p = pow(c, 2.0)/mu;

        // Calculate perturbation components
        W<T> Fr = thames::vector::geometry::dot3(F, er);
        W<T> Fh = thames::vector::geometry::dot3(F, eh);

        // Calculate non-dimensional quantities
        W<T> zeta = r/p;
        W<T> zetatilde = 1 + zeta;

        // Calculate time derivatives of the second and third elements
        W<T> p1dot = p2*((h - c)/pow(r, 2.0) - r/h*hwh*Fh) + 1.0/c*(r*drdt/c*p1 + zetatilde*p2 + zeta*cl)*(2.0*U - r*Fr) + r/mu*(zeta*p1 + zetatilde*sl)*edot;
        W<T> p2dot = p1*(r/h*hwh*Fh - (h-c)/pow(r, 2.0)) + 1.0/c*(r*drdt/c*p2 - zetatilde*p1 - zeta*sl)*(2.0*U - r*Fr) + r/mu*(zeta*p2 + zetatilde*cl)*edot;

        // Calculate time derivative of the generalised mean longitude
        W<T> Ldot = nu + (h - c)/pow(r, 2.0) - r/h*hwh*Fh + (r*drdt*c/pow(mu, 2.0)*zetatilde*alpha)*edot + 1.0/c*(1.0/alpha + alpha*(1.0 - r/a))*(2.0*U - r*Fr);

        // Calculate time derivatives of the remaining elements
        W<T> q1dot = r/(2.0*h)*Fh*(1.0 + pow(q1, 2.0) + pow(q2, 2.0))*sl;
        W<T> q2dot = r/(2.0*h)*Fh*(1.0 + pow(q1, 2.0) + pow(q2, 2.0))*cl;

        // Store derivatives
        geqoedot = {
            nudot,
            p1dot,
            p2dot,
            Ldot,
            q1dot,
            q2dot
        };

        // Return zero
        return 0;
    }

    template class GEqOEPropagatorPolynomialDynamics<double, taylor_polynomial>;
    template class GEqOEPropagatorPolynomialDynamics<double, chebyshev_polynomial>;

    template<class T, template<class> class P>
    GEqOEPropagatorPolynomial<T, P>::GEqOEPropagatorPolynomial(const T& mu, BasePerturbationPolynomial<T, P>* const perturbation, const DimensionalFactors<T>* factors) : BasePropagatorPolynomial<T, P>(mu, perturbation, factors, new GEqOEPropagatorPolynomialDynamics<T, P>(mu, perturbation, factors), GEQOE) {

    }

    template<class T, template<class> class P>
    GEqOEPropagatorPolynomial<T, P>::~GEqOEPropagatorPolynomial() {

    }

    template class GEqOEPropagatorPolynomial<double, taylor_polynomial>;
    template class GEqOEPropagatorPolynomial<double, chebyshev_polynomial>;

    #endif

}