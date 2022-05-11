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

#include "../../include/conversions/dimensional.h"
#include "../../include/propagators/cowell.h"
#include "../../include/propagators/options.h"
#include "../../include/perturbations/baseperturbation.h"
#include "../../include/vector/arithmeticoverloads.h"
#include "../../include/vector/geometry.h"

namespace thames::propagators {

    using thames::perturbations::baseperturbation::BasePerturbation;
    using namespace thames::vector::arithmeticoverloads;
    using thames::conversions::dimensional::DimensionalFactors;

    ///////////
    // Reals //
    ///////////

    template<class T>
    CowellPropagator<T>::CowellPropagator(const T& mu, const BasePerturbation<T>* perturbation, const DimensionalFactors<T>* factors) : BasePropagator<T>(factors), m_mu(mu/factors->grav), m_perturbation(perturbation) {

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
    std::array<T, 6> CowellPropagator<T>::propagate(T tstart, T tend, T tstep, std::array<T, 6> state, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype) const {
        // Check that input is Cartesian state
        if(statetype != thames::constants::statetypes::CARTESIAN)
            throw std::runtime_error("Unsupported state type");
        
        // Non-dimensionalise times
        tstart /= m_factors->time;
        tend /= m_factors->time;
        tstep /= m_factors->time;

        // Non-dimensionalise state
        state = thames::conversions::dimensional::cartesian_nondimensionalise(state, *m_factors);
        
        // Declare state derivative
        auto func = [this](const std::array<T, 6>& x, std::array<T, 6>& dxdt, const T t){return derivative(x, dxdt, t);};

        // Propagate according to the fixed flag
        if(options.fixed){
            // Declare stepper
            boost::numeric::odeint::runge_kutta4<std::array<T, 6>> stepper;

            // Propagate orbit
            boost::numeric::odeint::integrate_const(stepper, func, state, tstart, tend, tstep);
        } else {
            // Declare stepper
            boost::numeric::odeint::runge_kutta_cash_karp54<std::array<T, 6>> stepper;
            auto steppercontrolled = boost::numeric::odeint::make_controlled(options.atol, options.rtol, stepper);

            // Propagate orbit
            boost::numeric::odeint::integrate_adaptive(steppercontrolled, func, state, tstart, tend, tstep);
        }

        // Re-dimensionalise state
        state = thames::conversions::dimensional::cartesian_dimensionalise(state, *m_factors);

        // Return final state
        return state;
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
    std::vector<T> CowellPropagator<T>::propagate(T tstart, T tend, T tstep, std::vector<T> state, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype) const {
        // Check that input is Cartesian state
        if(statetype != thames::constants::statetypes::CARTESIAN)
            throw std::runtime_error("Unsupported state type");
        
        // Non-dimensionalise times
        tstart /= m_factors->time;
        tend /= m_factors->time;
        tstep /= m_factors->time;

        // Non-dimensionalise state
        state = thames::conversions::dimensional::cartesian_nondimensionalise(state, *m_factors);
        
        // Declare state derivative
        auto func = [this](const std::vector<T>& x, std::vector<T>& dxdt, const T t){return derivative(x, dxdt, t);};

        // Propagate according to the fixed flag
        if(options.fixed){
            // Declare stepper
            boost::numeric::odeint::runge_kutta4<std::vector<T>> stepper;

            // Propagate orbit
            boost::numeric::odeint::integrate_const(stepper, func, state, tstart, tend, tstep);
        } else {
            // Declare stepper
            boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<T>> stepper;
            auto steppercontrolled = boost::numeric::odeint::make_controlled(options.atol, options.rtol, stepper);

            // Propagate orbit
            boost::numeric::odeint::integrate_adaptive(steppercontrolled, func, state, tstart, tend, tstep);
        }

        // Re-dimensionalise state
        state = thames::conversions::dimensional::cartesian_dimensionalise(state, *m_factors);

        // Return final state
        return state;
    }

    template<class T>
    std::vector<std::vector<T>> CowellPropagator<T>::propagate(std::vector<T> tvector, T tstep, std::vector<T> state, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype) const {
        // Declare output vector
        std::vector<std::vector<T>> states(tvector.size());

        // Add initial state to output vector
        states[0] = state;

        // Propagate between times
        for(std::size_t ii=0; ii<tvector.size()-1; ii++){
            states[ii+1] = propagate(tvector[ii], tvector[ii+1], tstep, states[ii], options, statetype);
        }

        // Return states
        return states;
    }

    template class CowellPropagator<double>;

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::integrator;
    using namespace smartuq::polynomial;
    using thames::perturbations::baseperturbation::BasePerturbationPolynomial;

    template<class T, template<class> class P>
    CowellPropagatorPolynomialDynamics<T, P>::CowellPropagatorPolynomialDynamics(const T& mu, const BasePerturbationPolynomial<T, P>* perturbation) : smartuq::dynamics::base_dynamics<P<T>>("Cowell"), m_mu(mu), m_perturbation(perturbation) {

    }

    template<class T, template<class> class P>
    CowellPropagatorPolynomialDynamics<T, P>::~CowellPropagatorPolynomialDynamics() {

    }

    template<class T, template<class> class P>
    int CowellPropagatorPolynomialDynamics<T, P>::evaluate(const T& t, const std::vector<P<T>>& RV, std::vector<P<T>>& RVdot) const {
        // Extract Cartesian state vectors
        std::vector<P<T>> R = {RV[0], RV[1], RV[2]};
        std::vector<P<T>> V = {RV[3], RV[4], RV[5]};

        // Calculate range
        P<T> r = thames::vector::geometry::norm3(R);

        // Calculate perturbing acceleration
        std::vector<P<T>> F = m_perturbation->acceleration_total(t, R, V);

        // Calculate central body acceleration
        std::vector<P<T>> G = -m_mu/pow(r, 3)*R;

        // Calculate acceleration
        std::vector<P<T>> A = G + F;

        // Update derivative
        RVdot = {V[0], V[1], V[2], A[0], A[1], A[2]};

        // Return zero
        return 0;
    }

    template class CowellPropagatorPolynomialDynamics<double, taylor_polynomial>;
    template class CowellPropagatorPolynomialDynamics<double, chebyshev_polynomial>;

    template<class T, template<class> class P>
    CowellPropagatorPolynomial<T, P>::CowellPropagatorPolynomial(const T& mu, const BasePerturbationPolynomial<T, P>* perturbation, const DimensionalFactors<T>* factors) : BasePropagatorPolynomial<T, P>(factors), m_dyn(mu/factors->grav, perturbation) {

    }

    template<class T, template<class> class P>
    CowellPropagatorPolynomial<T, P>::~CowellPropagatorPolynomial() {
        
    }

    template<class T, template<class> class P>
    std::vector<P<T>> CowellPropagatorPolynomial<T, P>::propagate(T tstart, T tend, T tstep, std::vector<P<T>> state, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype) const {
        // Check that input is Cartesian state
        if(statetype != thames::constants::statetypes::CARTESIAN)
            throw std::runtime_error("Unsupported state type");
        
        // Non-dimensionalise times
        tstart /= m_factors->time;
        tend /= m_factors->time;
        tstep /= m_factors->time;

        // Non-dimensionalise state
        state = thames::conversions::dimensional::cartesian_nondimensionalise(state, *m_factors);
        
        // Calculate number of steps based on time step
        unsigned int nstep = (int) ceil((tend - tstart)/tstep);

        // Create final state vector
        std::vector<P<T>> statefinal(state);

        // Propagate according to the fixed flag
        if(options.fixed){
            // Create integrator
            rk4<P<T>> integrator(&m_dyn);

            // Integrate state
            integrator.integrate(tstart, tend, nstep, state, statefinal);  
        } else {
            // Create integrator
            rk45<P<T>> integrator(&m_dyn, options.atol, options.rtol);

            // Integrate state
            integrator.integrate(tstart, tend, nstep, state, statefinal);  
        }

        // Re-dimensionalise state
        statefinal = thames::conversions::dimensional::cartesian_dimensionalise(statefinal, *m_factors);

        // Return final state
        return statefinal;        
    }

    template<class T, template<class> class P>
    std::vector<std::vector<P<T>>> CowellPropagatorPolynomial<T, P>::propagate(std::vector<T> tvector, T tstep, std::vector<P<T>> state, thames::propagators::options::PropagatorOptions<T> options, thames::constants::statetypes::StateTypes statetype) const {
        // Declare output vector
        std::vector<std::vector<P<T>>> states(tvector.size());

        // Add initial state to output vector
        states[0] = state;

        // Propagate between times
        for(std::size_t ii=0; ii<tvector.size()-1; ii++){
            states[ii+1] = propagate(tvector[ii], tvector[ii+1], tstep, states[ii], options, statetype);
        }

        // Return states
        return states;
    }

    template class CowellPropagatorPolynomial<double, taylor_polynomial>;
    template class CowellPropagatorPolynomial<double, chebyshev_polynomial>;

    #endif

}