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
#include <vector>

#include <boost/numeric/odeint.hpp>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Integrators/rk4.h"
#include "../../external/smart-uq/include/Integrators/rk45.h"
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../../include/constants/statetypes.h"
#include "../../include/conversions/dimensional.h"
#include "../../include/conversions/universal.h"
#include "../../include/propagators/basepropagator.h"
#include "../../include/propagators/options.h"

namespace thames::propagators::basepropagator {

    using thames::constants::statetypes::StateTypes;
    using thames::constants::statetypes::CARTESIAN;
    using thames::propagators::options::PropagatorOptions;

    template<class T>
    BasePropagator<T>::BasePropagator(const T& mu, BasePerturbation<T>* const perturbation, const DimensionalFactors<T>* factors, const StateTypes propstatetype) : m_mu(mu), m_perturbation(perturbation), m_factors(factors), m_propstatetype(propstatetype) {

    }

    template<class T>
    BasePropagator<T>::~BasePropagator() {

    }

    ////////////
    // Arrays //
    ////////////

    template<class T>
    void BasePropagator<T>::derivative(const std::array<T, 6>& RV, std::array<T, 6>& RVdot, const T t) const {
        // Throw error if derivative is not implemented in dervied propagators
        throw std::runtime_error("Derivative must be defined");
    }

    template<class T>
    std::array<T, 6> BasePropagator<T>::propagate(T tstart, T tend, T tstep, std::array<T, 6> state, const PropagatorOptions<T> options, const StateTypes statetype) {
        // Check that input is Cartesian state
        if(statetype != CARTESIAN)
            throw std::runtime_error("Unsupported state type");

        // Set non-dimensional flag
        m_isNonDimensional = options.isNonDimensional;

        // Set perturbation non-dimensional flag
        m_perturbation->set_nondimensional(options.isNonDimensional);

        // Non-dimensionalise
        if (options.isNonDimensional) {
            tstart /= m_factors->time;
            tend /= m_factors->time;
            tstep /= m_factors->time;
            state = thames::conversions::dimensional::cartesian_nondimensionalise(state, *m_factors);
        }

        // Calculate gravitational parameters
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;

        // Convert state
        state = thames::conversions::universal::convert_state(tstart, state, mu, statetype, m_propstatetype, m_perturbation);

        // Declare state derivative
        auto func = [this](const std::array<T, 6>& x, std::array<T, 6>& dxdt, const T t){return derivative(x, dxdt, t);};

        // Propagate according to the fixed flag
        if(options.isfixedStep){
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

        // Convert state
        state = thames::conversions::universal::convert_state(tend, state, mu, m_propstatetype, statetype, m_perturbation);

        // Re-dimensionalise
        if (options.isNonDimensional) {
            state = thames::conversions::dimensional::cartesian_dimensionalise(state, *m_factors);
        }

        // Return final state
        return state;

    }

    /////////////
    // Vectors //
    /////////////

    template<class T>
    void BasePropagator<T>::derivative(const std::vector<T>& RV, std::vector<T>& RVdot, const T t) const {
        // Throw error if derivative is not implemented in dervied propagators
        throw std::runtime_error("Derivative must be defined");
    }

    template<class T>
    std::vector<T> BasePropagator<T>::propagate(T tstart, T tend, T tstep, std::vector<T> state, const PropagatorOptions<T> options, const StateTypes statetype) {
        // Check that input is Cartesian state
        if(statetype != CARTESIAN)
            throw std::runtime_error("Unsupported state type");

        // Set non-dimensional flag
        m_isNonDimensional = options.isNonDimensional;

        // Set perturbation non-dimensional flag
        m_perturbation->set_nondimensional(options.isNonDimensional);

        // Non-dimensionalise
        if (options.isNonDimensional) {
            tstart /= m_factors->time;
            tend /= m_factors->time;
            tstep /= m_factors->time;
            state = thames::conversions::dimensional::cartesian_nondimensionalise(state, *m_factors);
        }

        // Calculate gravitational parameters
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;

        // Convert state
        state = thames::conversions::universal::convert_state(tstart, state, mu, statetype, m_propstatetype, m_perturbation);

        // Declare state derivative
        auto func = [this](const std::vector<T>& x, std::vector<T>& dxdt, const T t){return derivative(x, dxdt, t);};

        // Propagate according to the fixed flag
        if(options.isfixedStep){
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

        // Convert state
        state = thames::conversions::universal::convert_state(tend, state, mu, m_propstatetype, statetype, m_perturbation);

        // Re-dimensionalise
        if (options.isNonDimensional) {
            state = thames::conversions::dimensional::cartesian_dimensionalise(state, *m_factors);
        }

        // Return final state
        return state;

    }

    template class BasePropagator<double>;

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::integrator;
    using namespace smartuq::polynomial;

    template<class T, template<class> class P>
    BasePropagatorPolynomialDynamics<T, P>::BasePropagatorPolynomialDynamics(std::string name, const T& mu, BasePerturbationPolynomial<T, P>* const perturbation, const DimensionalFactors<T>* factors) : smartuq::dynamics::base_dynamics<P<T>>(name), m_mu(mu), m_perturbation(perturbation), m_factors(factors) {

    }

    template<class T, template<class> class P>
    BasePropagatorPolynomialDynamics<T, P>::~BasePropagatorPolynomialDynamics() {

    }

    template<class T, template<class> class P>
    int BasePropagatorPolynomialDynamics<T, P>::evaluate(const T& t, const std::vector<P<T>>& x, std::vector<P<T>>& dxdt) const {
        // Throw error if derivative is not implemented in dervied propagators
        throw std::runtime_error("Derivative must be defined");
    }

    template class BasePropagatorPolynomialDynamics<double, taylor_polynomial>;
    template class BasePropagatorPolynomialDynamics<double, chebyshev_polynomial>;

    template<class T, template<class> class P>
    BasePropagatorPolynomial<T, P>::BasePropagatorPolynomial(const T& mu, BasePerturbationPolynomial<T, P>* const perturbation, const DimensionalFactors<T>* factors, BasePropagatorPolynomialDynamics<T, P>* const dyn, const StateTypes propstatetype) : m_mu(mu), m_perturbation(perturbation), m_factors(factors), m_dyn(dyn), m_propstatetype(propstatetype) {

    }

    template<class T, template<class> class P>
    std::vector<P<T>> BasePropagatorPolynomial<T, P>::propagate(T tstart, T tend, T tstep, std::vector<P<T>> state, PropagatorOptions<T> options, StateTypes statetype) {
        // Check that input is Cartesian state
        if(statetype != thames::constants::statetypes::CARTESIAN)
            throw std::runtime_error("Unsupported state type");

        // Set non-dimensional flag
        m_dyn->m_isNonDimensional = options.isNonDimensional;

        // Set perturbation non-dimensional flag
        m_perturbation->set_nondimensional(options.isNonDimensional);
        
        // Non-dimensionalise
        if (options.isNonDimensional) {
            tstart /= m_factors->time;
            tend /= m_factors->time;
            tstep /= m_factors->time;
            state = thames::conversions::dimensional::cartesian_nondimensionalise(state, *m_factors);
        }

        // Calculate gravitational parameters
        const T mu = (options.isNonDimensional) ? m_mu/m_factors->grav : m_mu;

        // Convert state
        state = thames::conversions::universal::convert_state(tstart, state, mu, statetype, m_propstatetype, m_perturbation);
        
        // Calculate number of steps based on time step
        unsigned int nstep = (int) ceil((tend - tstart)/tstep);

        // Create final state vector
        std::vector<P<T>> statefinal(state);

        // Propagate according to the fixed flag
        if(options.isfixedStep){
            // Create integrator
            rk4<P<T>> integrator(m_dyn);

            // Integrate state
            integrator.integrate(tstart, tend, nstep, state, statefinal);  
        } else {
            // Create integrator
            rk45<P<T>> integrator(m_dyn, options.atol, options.rtol);

            // Integrate state
            integrator.integrate(tstart, tend, nstep, state, statefinal);  
        }

        // Convert state
        statefinal = thames::conversions::universal::convert_state(tend, statefinal, mu, m_propstatetype, statetype, m_perturbation);
        
        // Re-dimensionalise
        if (options.isNonDimensional) {
            statefinal = thames::conversions::dimensional::cartesian_dimensionalise(statefinal, *m_factors);
        }

        // Return final state
        return statefinal;        
    }

    template class BasePropagatorPolynomial<double, taylor_polynomial>;
    template class BasePropagatorPolynomial<double, chebyshev_polynomial>;

    #endif

}