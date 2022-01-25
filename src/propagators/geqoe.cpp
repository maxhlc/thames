#include <array>
#include <cmath>
#include <vector>

#include <boost/numeric/odeint.hpp>

#include "geqoe.h"
#include "../conversions/geqoe.h"
#include "../perturbations/baseperturbation.h"
#include "../vector/arithmeticoverloads.h"
#include "../vector/geometry.h"

using namespace thames::perturbations::baseperturbation;
using namespace thames::vector::arithmeticoverloads;

namespace thames::propagators::geqoe{

    template<class T>
    GEqOEPropagator<T>::GEqOEPropagator(const T& mu, BasePerturbation<T>* perturbation){
        // Store gravitational parameter
        m_mu = mu;

        // Store perturbation object
        m_perturbation = perturbation;
    }

    ////////////
    // Arrays //
    ////////////

    template<class T>
    void GEqOEPropagator<T>::derivative(const std::array<T, 6>& geqoe, std::array<T, 6>& geqoedot, const T t) const {
        // Extract elements
        T nu = geqoe[0];
        T p1 = geqoe[1];
        T p2 = geqoe[2];
        T q1 = geqoe[4];
        T q2 = geqoe[5];

        // Calculate Cartesian state and extract vectors
        std::array<T, 6> RV = thames::conversions::geqoe::geqoe_to_cartesian(t, geqoe, m_mu, m_perturbation);
        std::array<T, 3> R = {RV[0], RV[1], RV[2]};
        std::array<T, 3> V = {RV[3], RV[4], RV[5]};

        // Calculate range and range rate
        T r = thames::vector::geometry::norm3(R);
        T drdt = thames::vector::geometry::dot3(R, V)/r;

        // Calculate perturbations
        T U = m_perturbation->potential(t, R);
        T Ut = m_perturbation->potential_derivative(t, R, V);
        std::array<T, 3> F = m_perturbation->acceleration_total(t, R, V);
        std::array<T, 3> P = m_perturbation->acceleration_nonpotential(t, R, V);

        // Calculate time derivative of total energy
        T edot = Ut + thames::vector::geometry::dot3(P, V);

        // Calculate time derivative of nu
        T nudot = -3.0*pow(nu/pow(m_mu, 2.0), 1.0/3.0)*edot;

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

        // Calculate equinoctial reference frame velocity components
        T hwh = q1*cl - q2*sl;

        // Calculate angular momentum
        std::array<T, 3> H = thames::vector::geometry::cross3(R, V);
        T h = thames::vector::geometry::norm3(H);
        std::array<T, 3> eh = H/h;

        // Calculate the effective potential energy
        T ueff = pow(h, 2.0)/(2.0*pow(r, 2.0)) + U;

        // Calculate the generalised angular momentum
        T c = sqrt(2.0*pow(r, 2.0)*ueff);

        // Calculate the generalised semi-latus rectum
        T p = pow(c, 2.0)/m_mu;

        // Calculate perturbation components
        T Fr = thames::vector::geometry::dot3(F, er);
        T Fh = thames::vector::geometry::dot3(F, eh);

        // Calculate non-dimensional quantities
        T zeta = r/p;
        T zetatilde = 1 + zeta;

        // Calculate time derivatives of the second and third elements
        T p1dot = p2*((h - c)/pow(r, 2.0) - r/h*hwh*Fh) + 1.0/c*(r*drdt/c*p1 + zetatilde*p2 + zeta*cl)*(2.0*U - r*Fr) + r/m_mu*(zeta*p1 + zetatilde*sl)*edot;
        T p2dot = p1*(r/h*hwh*Fh - (h-c)/pow(r, 2.0)) + 1.0/c*(r*drdt/c*p2 - zetatilde*p1 - zeta*sl)*(2.0*U - r*Fr) + r/m_mu*(zeta*p2 + zetatilde*cl)*edot;

        // Calculate generalised semi-major axis
        T a = pow(m_mu/pow(nu, 2.0), 1.0/3.0);

        // Calculate alpha
        T alpha = 1.0/(1.0 + sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0)));

        // Calculate time derivative of the generalised mean longitude
        T Ldot = nu + (h - c)/pow(r, 2.0) - r/h*hwh*Fh + (r*drdt*c/pow(m_mu, 2.0)*zetatilde*alpha)*edot + 1.0/c*(1.0/alpha + alpha*(1.0 - r/a))*(2.0*U - r*Fr);

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

    template<class T>
    std::array<T, 6> GEqOEPropagator<T>::propagate(T tstart, T tend, T tstep, std::array<T, 6> RV, T atol, T rtol) const {
        // Transform initial state
        std::array<T, 6> geqoe = thames::conversions::geqoe::cartesian_to_geqoe(tstart, RV, m_mu, m_perturbation);

        // Declare stepper
        boost::numeric::odeint::runge_kutta_cash_karp54<std::array<T, 6>> stepper;
        auto steppercontrolled = boost::numeric::odeint::make_controlled(atol, rtol, stepper);

        // Propagate orbit
        boost::numeric::odeint::integrate_adaptive(steppercontrolled, [this](const std::array<T, 6>& x, std::array<T, 6>& dxdt, const T t){return derivative(x, dxdt, t);}, geqoe, tstart, tend, tstep);

        // Transform final state
        RV = thames::conversions::geqoe::geqoe_to_cartesian<T>(tend, geqoe, m_mu, m_perturbation);

        // Return final state
        return RV;
    }

    /////////////
    // Vectors //
    /////////////

    template<class T>
    void GEqOEPropagator<T>::derivative(const std::vector<T>& geqoe, std::vector<T>& geqoedot, const T t) const {
        // Extract elements
        T nu = geqoe[0];
        T p1 = geqoe[1];
        T p2 = geqoe[2];
        T q1 = geqoe[4];
        T q2 = geqoe[5];

        // Calculate Cartesian state and extract vectors
        std::vector<T> RV = thames::conversions::geqoe::geqoe_to_cartesian(t, geqoe, m_mu, m_perturbation);
        std::vector<T> R = {RV[0], RV[1], RV[2]};
        std::vector<T> V = {RV[3], RV[4], RV[5]};

        // Calculate range and range rate
        T r = thames::vector::geometry::norm3(R);
        T drdt = thames::vector::geometry::dot3(R, V)/r;

        // Calculate perturbations
        T U = m_perturbation->potential(t, R);
        T Ut = m_perturbation->potential_derivative(t, R, V);
        std::vector<T> F = m_perturbation->acceleration_total(t, R, V);
        std::vector<T> P = m_perturbation->acceleration_nonpotential(t, R, V);

        // Calculate time derivative of total energy
        T edot = Ut + thames::vector::geometry::dot3(P, V);

        // Calculate time derivative of nu
        T nudot = -3.0*pow(nu/pow(m_mu, 2.0), 1.0/3.0)*edot;

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

        // Calculate equinoctial reference frame velocity components
        T hwh = q1*cl - q2*sl;

        // Calculate angular momentum
        std::vector<T> H = thames::vector::geometry::cross3(R, V);
        T h = thames::vector::geometry::norm3(H);
        std::vector<T> eh = H/h;

        // Calculate the effective potential energy
        T ueff = pow(h, 2.0)/(2.0*pow(r, 2.0)) + U;

        // Calculate the generalised angular momentum
        T c = sqrt(2.0*pow(r, 2.0)*ueff);

        // Calculate the generalised semi-latus rectum
        T p = pow(c, 2.0)/m_mu;

        // Calculate perturbation components
        T Fr = thames::vector::geometry::dot3(F, er);
        T Fh = thames::vector::geometry::dot3(F, eh);

        // Calculate non-dimensional quantities
        T zeta = r/p;
        T zetatilde = 1 + zeta;

        // Calculate time derivatives of the second and third elements
        T p1dot = p2*((h - c)/pow(r, 2.0) - r/h*hwh*Fh) + 1.0/c*(r*drdt/c*p1 + zetatilde*p2 + zeta*cl)*(2.0*U - r*Fr) + r/m_mu*(zeta*p1 + zetatilde*sl)*edot;
        T p2dot = p1*(r/h*hwh*Fh - (h-c)/pow(r, 2.0)) + 1.0/c*(r*drdt/c*p2 - zetatilde*p1 - zeta*sl)*(2.0*U - r*Fr) + r/m_mu*(zeta*p2 + zetatilde*cl)*edot;

        // Calculate generalised semi-major axis
        T a = pow(m_mu/pow(nu, 2.0), 1.0/3.0);

        // Calculate alpha
        T alpha = 1.0/(1.0 + sqrt(1.0 - pow(p1, 2.0) - pow(p2, 2.0)));

        // Calculate time derivative of the generalised mean longitude
        T Ldot = nu + (h - c)/pow(r, 2.0) - r/h*hwh*Fh + (r*drdt*c/pow(m_mu, 2.0)*zetatilde*alpha)*edot + 1.0/c*(1.0/alpha + alpha*(1.0 - r/a))*(2.0*U - r*Fr);

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

    template<class T>
    std::vector<T> GEqOEPropagator<T>::propagate(T tstart, T tend, T tstep, std::vector<T> RV, T atol, T rtol) const {
        // Transform initial state
        std::vector<T> geqoe = thames::conversions::geqoe::cartesian_to_geqoe<T>(tstart, RV, m_mu, m_perturbation);

        // Declare stepper
        boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<T>> stepper;
        auto steppercontrolled = boost::numeric::odeint::make_controlled(atol, rtol, stepper);

        // Propagate orbit
        boost::numeric::odeint::integrate_adaptive(steppercontrolled, [this](const std::vector<T>& x, std::vector<T>& dxdt, const T t){return derivative(x, dxdt, t);}, geqoe, tstart, tend, tstep);

        // Transform final state
        RV = thames::conversions::geqoe::geqoe_to_cartesian<T>(tend, geqoe, m_mu, m_perturbation);

        // Return final state
        return RV;
    }

    template class GEqOEPropagator<double>;

}