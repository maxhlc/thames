#include <array>
#include <cmath>
#include <vector>

#include <boost/numeric/odeint.hpp>

#include "cowell.h"
#include "../perturbations/baseperturbation.h"
#include "../vector/arithmeticoverloads.h"
#include "../vector/geometry.h"

using namespace thames::perturbations::baseperturbation;
using namespace thames::vector::arithmeticoverloads;

namespace thames::propagators {

    template<class T>
    CowellPropagator<T>::CowellPropagator(const T& mu, const BasePerturbation<T>* perturbation) : m_mu(mu), m_perturbation(perturbation) {

    }

    ////////////
    // Arrays //
    ////////////

    template<class T>
    void CowellPropagator<T>::derivative(const std::array<T, 6>& RV, std::array<T, 6>& RVdot, const T t) const {
        // Extract Cartesian state vectors
        std::array<T, 3> R = {RV[0], RV[1], RV[2]};
        std::array<T, 3> V = {RV[3], RV[4], RV[5]};

        // Calculate range
        T r = thames::vector::geometry::norm3(R);

        // Calculate perturbing acceleration
        std::array<T, 3> F = m_perturbation->acceleration_total(t, R, V);

        // Calculate central body acceleration
        std::array<T, 3> G = -m_mu/pow(r, 3.0)*R;

        // Calculate acceleration
        std::array<T, 3> A = G + F;

        // Store state derivative
        for(unsigned int ii=0; ii<3; ii++){
            RVdot[ii] = V[ii];
            RVdot[ii+3] = A[ii];
        }
    }

    template<class T>
    std::array<T, 6> CowellPropagator<T>::propagate(T tstart, T tend, T tstep, std::array<T, 6> RV,  T atol, T rtol) const {
        // Declare stepper
        boost::numeric::odeint::runge_kutta_cash_karp54<std::array<T, 6>> stepper;
        auto steppercontrolled = boost::numeric::odeint::make_controlled(atol, rtol, stepper);

        // Propagate orbit
        boost::numeric::odeint::integrate_adaptive(steppercontrolled, [this](const std::array<T, 6>& x, std::array<T, 6>& dxdt, const T t){return derivative(x, dxdt, t);}, RV, tstart, tend, tstep);

        // Return final state
        return RV;
    }

    /////////////
    // Vectors //
    /////////////

    template<class T>
    void CowellPropagator<T>::derivative(const std::vector<T>& RV, std::vector<T>& RVdot, const T t) const {
        // Extract Cartesian state vectors
        std::vector<T> R = {RV[0], RV[1], RV[2]};
        std::vector<T> V = {RV[3], RV[4], RV[5]};

        // Calculate range
        T r = thames::vector::geometry::norm3(R);

        // Calculate perturbing acceleration
        std::vector<T> F = m_perturbation->acceleration_total(t, R, V);

        // Calculate central body acceleration
        std::vector<T> G = -m_mu/pow(r, 3.0)*R;

        // Calculate acceleration
        std::vector<T> A = G + F;

        // Store state derivative
        for(unsigned int ii=0; ii<3; ii++){
            RVdot[ii] = V[ii];
            RVdot[ii+3] = A[ii];
        }
    }

    template<class T>
    std::vector<T> CowellPropagator<T>::propagate(T tstart, T tend, T tstep, std::vector<T> RV, T atol, T rtol) const {
        // Declare stepper
        boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<T>> stepper;
        auto steppercontrolled = boost::numeric::odeint::make_controlled(atol, rtol, stepper);

        // Propagate orbit
        boost::numeric::odeint::integrate_adaptive(steppercontrolled, [this](const std::vector<T>& x, std::vector<T>& dxdt, const T t){return derivative(x, dxdt, t);}, RV, tstart, tend, tstep);

        // Return final state
        return RV;
    }

    template class CowellPropagator<double>;

}