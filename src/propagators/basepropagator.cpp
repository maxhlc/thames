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
#include <memory>
#include <vector>

#include <boost/numeric/odeint.hpp>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Integrators/rk4.h"
#include "../../external/smart-uq/include/Integrators/rk45.h"
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../../include/constants/statetypes.h"
#include "../../include/conversions/cartesian.h"
#include "../../include/conversions/dimensional.h"
#include "../../include/conversions/universal.h"
#include "../../include/propagators/basepropagator.h"
#include "../../include/settings/settings.h"
#include "../../include/util/polynomials.h"

namespace thames::propagators::basepropagator {

    using thames::constants::statetypes::StateTypes;
    using thames::constants::statetypes::CARTESIAN;
    using thames::settings::PropagatorParameters;

    template<class T>
    BasePropagator<T>::BasePropagator(const T& mu, const std::shared_ptr<BasePerturbation<T>> perturbation, const std::shared_ptr<DimensionalFactors<T>> factors, const StateTypes propstatetype) : m_mu(mu), m_perturbation(perturbation), m_factors(factors), m_propstatetype(propstatetype) {

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
    std::array<T, 6> BasePropagator<T>::propagate(T tstart, T tend, T tstep, std::array<T, 6> state, const PropagatorParameters<T> options, const StateTypes statetype) {
        // Check that input is Cartesian state
        if(statetype != CARTESIAN)
            throw std::runtime_error("Unsupported state type");

        // Set non-dimensional flag
        m_isNonDimensional = options.isNonDimensional;

        // Set perturbation non-dimensional flag
        m_perturbation->set_nondimensional(options.isNonDimensional);

        // Non-dimensionalise
        if (options.isNonDimensional) {
            // Update factors
            *m_factors = thames::conversions::dimensional::calculate_factors(state, m_mu);

            // Scale times
            tstart /= m_factors->time;
            tend /= m_factors->time;
            tstep /= m_factors->time;

            // Scale state
            state = thames::conversions::dimensional::cartesian_nondimensionalise(state, *m_factors);
        }

        // Calculate gravitational parameters
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;

        // Convert state
        state = thames::conversions::universal::convert_state<T>(tstart, state, mu, statetype, m_propstatetype, m_perturbation);

        // Declare state derivative
        auto func = [this](const std::array<T, 6>& x, std::array<T, 6>& dxdt, const T t){return derivative(x, dxdt, t);};

        // Propagate according to the fixed flag
        if(options.isFixedStep){
            // Declare stepper
            boost::numeric::odeint::runge_kutta4<std::array<T, 6>> stepper;

            // Propagate orbit
            boost::numeric::odeint::integrate_const(stepper, func, state, tstart, tend, tstep);
        } else {
            // Declare stepper
            boost::numeric::odeint::runge_kutta_cash_karp54<std::array<T, 6>> stepper;
            auto steppercontrolled = boost::numeric::odeint::make_controlled(options.absoluteTolerance, options.relativeTolerance, stepper);

            // Propagate orbit
            boost::numeric::odeint::integrate_adaptive(steppercontrolled, func, state, tstart, tend, tstep);
        }

        // Convert state
        state = thames::conversions::universal::convert_state<T>(tend, state, mu, m_propstatetype, statetype, m_perturbation);

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
    std::vector<T> BasePropagator<T>::propagate(T tstart, T tend, T tstep, std::vector<T> state, const PropagatorParameters<T> options, const StateTypes statetype) {
        // Check that input is Cartesian state
        if(statetype != CARTESIAN)
            throw std::runtime_error("Unsupported state type");

        // Set non-dimensional flag
        m_isNonDimensional = options.isNonDimensional;

        // Set perturbation non-dimensional flag
        m_perturbation->set_nondimensional(options.isNonDimensional);

        // Non-dimensionalise
        if (options.isNonDimensional) {
            // Update factors
            *m_factors = thames::conversions::dimensional::calculate_factors(state, m_mu);

            // Scale times
            tstart /= m_factors->time;
            tend /= m_factors->time;
            tstep /= m_factors->time;

            // Scale state
            state = thames::conversions::dimensional::cartesian_nondimensionalise(state, *m_factors);
        }

        // Calculate gravitational parameters
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;

        // Convert state
        state = thames::conversions::universal::convert_state<T>(tstart, state, mu, statetype, m_propstatetype, m_perturbation);

        // Declare state derivative
        auto func = [this](const std::vector<T>& x, std::vector<T>& dxdt, const T t){return derivative(x, dxdt, t);};

        // Propagate according to the fixed flag
        if(options.isFixedStep){
            // Declare stepper
            boost::numeric::odeint::runge_kutta4<std::vector<T>> stepper;

            // Propagate orbit
            boost::numeric::odeint::integrate_const(stepper, func, state, tstart, tend, tstep);
        } else {
            // Declare stepper
            boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<T>> stepper;
            auto steppercontrolled = boost::numeric::odeint::make_controlled(options.absoluteTolerance, options.relativeTolerance, stepper);

            // Propagate orbit
            boost::numeric::odeint::integrate_adaptive(steppercontrolled, func, state, tstart, tend, tstep);
        }

        // Convert state
        state = thames::conversions::universal::convert_state<T>(tend, state, mu, m_propstatetype, statetype, m_perturbation);

        // Re-dimensionalise
        if (options.isNonDimensional) {
            state = thames::conversions::dimensional::cartesian_dimensionalise(state, *m_factors);
        }

        // Return final state
        return state;

    }

    template<class T>
    std::vector<std::vector<T>> BasePropagator<T>::propagate(const std::vector<T> tvec, const T tstep, const std::vector<T> state, const PropagatorParameters<T> options, const StateTypes statetype) {
        // Declare output vectors
        std::vector<std::vector<T>> states_propagated(tvec.size());

        // Append initial state to output
        states_propagated[0] = state;

        // Propagate state between times
        for (std::size_t ii = 0; ii < tvec.size() - 1; ii++) {
            states_propagated[ii+1] = propagate(tvec[ii], tvec[ii+1], tstep, states_propagated[ii], options, statetype);
        }

        // Return output vector
        return states_propagated;
    }

    template<class T>
    std::vector<std::vector<T>> BasePropagator<T>::propagate(const T tstart, const T tend, const T tstep, const std::vector<std::vector<T>> states, const PropagatorParameters<T> options, const StateTypes statetype) {
        // Declare output states
        std::vector<std::vector<T>> states_propagated(states);

        // Iterate through states
        for (std::vector<T>& state : states_propagated) {
            state = propagate(tstart, tend, tstep, state, options, statetype);
        }

        // Return states
        return states_propagated;
    }

    template<class T>
    std::vector<std::vector<std::vector<T>>> BasePropagator<T>::propagate(const std::vector<T> tvec, const T tstep, const std::vector<std::vector<T>> states, const PropagatorParameters<T> options, const StateTypes statetype) {
        // Declare output vectors
        std::vector<std::vector<std::vector<T>>> states_propagated(tvec.size());

        // Append initial state to output
        states_propagated[0] = states;

        // Propagate state between times
        for (std::size_t ii = 0; ii < tvec.size() - 1; ii++) {
            states_propagated[ii+1] = propagate(tvec[ii], tvec[ii+1], tstep, states_propagated[ii], options, statetype);
        }

        // Return output vector
        return states_propagated;
    }

    template class BasePropagator<double>;

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::integrator;
    using namespace smartuq::polynomial;

    template<class T, template<class> class P>
    BasePropagatorPolynomialDynamics<T, P>::BasePropagatorPolynomialDynamics(std::string name, const T& mu, const std::shared_ptr<BasePerturbationPolynomial<T, P>> perturbation, const std::shared_ptr<const DimensionalFactors<T>> factors) : smartuq::dynamics::base_dynamics<P<T>>(name), m_mu(mu), m_perturbation(perturbation), m_factors(factors) {

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
    BasePropagatorPolynomial<T, P>::BasePropagatorPolynomial(const T& mu, const std::shared_ptr<BasePerturbationPolynomial<T, P>> perturbation, const std::shared_ptr<DimensionalFactors<T>> factors, const std::shared_ptr<BasePropagatorPolynomialDynamics<T, P>> dyn, const StateTypes propstatetype) : m_mu(mu), m_perturbation(perturbation), m_factors(factors), m_dyn(dyn), m_propstatetype(propstatetype) {

    }

    template<class T, template<class> class P>
    std::vector<P<T>> BasePropagatorPolynomial<T, P>::propagate(T tstart, T tend, T tstep, std::vector<P<T>> state, const PropagatorParameters<T> options, const StateTypes statetype) {
        // Check that input is Cartesian state
        if(statetype != thames::constants::statetypes::CARTESIAN)
            throw std::runtime_error("Unsupported state type");

        // Set non-dimensional flag
        m_dyn->m_isNonDimensional = options.isNonDimensional;

        // Set perturbation non-dimensional flag
        m_perturbation->set_nondimensional(options.isNonDimensional);
        
        // Non-dimensionalise
        if (options.isNonDimensional) {
            // Update factors
            *m_factors = thames::conversions::dimensional::calculate_factors(state, m_mu);

            // Scale times
            tstart /= m_factors->time;
            tend /= m_factors->time;
            tstep /= m_factors->time;

            // Scale state
            state = thames::conversions::dimensional::cartesian_nondimensionalise(state, *m_factors);
        }

        // Calculate gravitational parameters
        const T mu = (options.isNonDimensional) ? m_mu/m_factors->grav : m_mu;

        // Convert state
        state = thames::conversions::universal::convert_state<T, P>(tstart, state, mu, statetype, m_propstatetype, m_perturbation);
        
        // Calculate number of steps based on time step
        unsigned int nstep = (int) ceil((tend - tstart)/tstep);

        // Create final state vector
        std::vector<P<T>> statefinal(state);

        // Propagate according to the fixed flag
        if(options.isFixedStep){
            // Create integrator
            rk4<P<T>> integrator(m_dyn.get());

            // Integrate state
            integrator.integrate(tstart, tend, nstep, state, statefinal);  
        } else {
            // Create integrator
            rk45<P<T>> integrator(m_dyn.get(), options.absoluteTolerance, options.relativeTolerance);

            // Integrate state
            integrator.integrate(tstart, tend, nstep, state, statefinal);  
        }

        // Convert state
        statefinal = thames::conversions::universal::convert_state<T, P>(tend, statefinal, mu, m_propstatetype, statetype, m_perturbation);
        
        // Re-dimensionalise
        if (options.isNonDimensional) {
            statefinal = thames::conversions::dimensional::cartesian_dimensionalise(statefinal, *m_factors);
        }

        // Return final state
        return statefinal;        
    }

    template<class T, template <class> class P>
    std::vector<std::vector<P<T>>> BasePropagatorPolynomial<T, P>::propagate(const std::vector<T> tvec, const T tstep, const std::vector<P<T>> state, const PropagatorParameters<T> options, const StateTypes statetype) {
        // Declare output vectors
        std::vector<std::vector<P<T>>> states_propagated(tvec.size());

        // Append initial states to output
        states_propagated[0] = state;

        // Propagate state between times
        for (std::size_t ii = 0; ii < tvec.size() - 1; ii++) {
            states_propagated[ii+1] = propagate(tvec[ii], tvec[ii+1], tstep, states_propagated[ii], options, statetype);
        }

        // Return output vector
        return states_propagated;
    }

    template<class T, template <class> class P>
    std::vector<std::vector<T>> BasePropagatorPolynomial<T, P>::propagate(const T tstart, const T tend, const T tstep, std::vector<std::vector<T>> states, const PropagatorParameters<T> options, const StateTypes statetype, const unsigned int degree) {
        // Check that input is Cartesian state
        if(statetype != thames::constants::statetypes::CARTESIAN)
            throw std::runtime_error("Unsupported state type");

        // Generate polynomials
        std::vector<P<T>> statepolynomial;
        std::vector<T> lower, upper;
        thames::conversions::cartesian::cartesian_to_polynomial(states, degree, statepolynomial, lower, upper);

        // Calculate sample points
        std::vector<std::vector<T>> samples = thames::conversions::cartesian::state_to_sample(states, lower, upper);

        // Propagate polynomials
        statepolynomial = propagate(tstart, tend, tstep, statepolynomial, options, statetype);
 
        // Sample polynomials
        states = thames::util::polynomials::evaluate_polynomials(statepolynomial, samples);

        // Return propgated states
        return states;
    }

    template<class T, template <class> class P>
    std::vector<std::vector<std::vector<T>>> BasePropagatorPolynomial<T, P>::propagate(const std::vector<T> tvec, const T tstep, const std::vector<std::vector<T>> states, const PropagatorParameters<T> options, const StateTypes statetype, const unsigned int degree) {
        // Check that input is Cartesian state
        if(statetype != thames::constants::statetypes::CARTESIAN)
            throw std::runtime_error("Unsupported state type");

        // Declare output vectors
        std::vector<std::vector<std::vector<T>>> states_propagated(tvec.size());

        // Append initial states to output
        states_propagated[0] = states;

        // Generate polynomials
        std::vector<P<T>> statepolynomial;
        std::vector<T> lower, upper;
        thames::conversions::cartesian::cartesian_to_polynomial(states, degree, statepolynomial, lower, upper);

        // Calculate sample points
        std::vector<std::vector<T>> samples = thames::conversions::cartesian::state_to_sample(states, lower, upper);

        // Propagate state between times
        for (std::size_t ii = 0; ii < tvec.size() - 1; ii++) {
            // Update polynomials
            statepolynomial = propagate(tvec[ii], tvec[ii+1], tstep, statepolynomial, options, statetype);

            // Sample polynomials and store
            states_propagated[ii+1] = thames::util::polynomials::evaluate_polynomials(statepolynomial, samples);
        }

        // Return propagated states
        return states_propagated;        
    }

    template class BasePropagatorPolynomial<double, taylor_polynomial>;
    template class BasePropagatorPolynomial<double, chebyshev_polynomial>;

    #endif

}