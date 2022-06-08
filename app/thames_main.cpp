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

#include <chrono>
#include <stdexcept>

#include "../include/thames.h"

template<class T>
thames::settings::Parameters<T> propagate(const thames::settings::Parameters<T>& parameters) {
    // Load constants
    T J2 = thames::constants::earth::J2;
    T mu = thames::constants::earth::mu;
    T radius = thames::constants::earth::radius;
    T w = thames::constants::earth::w;

    // Declare factors
    auto factors = std::make_shared<thames::conversions::dimensional::DimensionalFactors<T>>();

    // Set up perturbations
    auto perturbation = std::make_shared<thames::perturbations::perturbationcombiner::PerturbationCombiner<T>>(factors);

    // Set up atmosphere model
    if (parameters.perturbation.atmosphere.isEnabled) {
        // Select atmosphere model
        if (parameters.perturbation.atmosphere.model == "USSA76") {
            auto atmospheremodel = thames::perturbations::atmosphere::models::USSA76;
            auto atmosphereperturbation = std::make_shared<thames::perturbations::atmosphere::drag::Drag<T>>(radius, w, parameters.spacecraft.Cd, parameters.spacecraft.dragArea, parameters.spacecraft.mass, atmospheremodel, factors);
            perturbation->add_model(atmosphereperturbation);
        } else if (parameters.perturbation.atmosphere.model == "Wertz") {
            auto atmospheremodel = thames::perturbations::atmosphere::models::WERTZ;
            auto atmosphereperturbation = std::make_shared<thames::perturbations::atmosphere::drag::Drag<T>>(radius, w, parameters.spacecraft.Cd, parameters.spacecraft.dragArea, parameters.spacecraft.mass, atmospheremodel, factors);
            perturbation->add_model(atmosphereperturbation);
        } else {
            throw std::runtime_error("Unsupported atmosphere model requested");
        }        
    }

    // Set up geopotential model
    if (parameters.perturbation.geopotential.isEnabled) {
        // Select geopotential model
        if (parameters.perturbation.geopotential.model == "J2") {
            auto geopotentialmodel = std::make_shared<thames::perturbations::geopotential::J2<T>>(mu, J2, radius, factors);
            perturbation->add_model(geopotentialmodel);
        } else {
            throw std::runtime_error("Unsupported geopotential model requested");
        }
    }

    // Import states
    std::vector<std::vector<T>> states = parameters.states[0].states;
    std::vector<std::vector<T>> states_propagated;
    T tstart = parameters.propagator.startTime;
    T tend = parameters.propagator.endTime;
    T tstep = parameters.propagator.timeStep;

    // Import state type
    thames::constants::statetypes::StateTypes statetype;
    if (parameters.states[0].statetype == "Cartesian") {
        statetype = thames::constants::statetypes::CARTESIAN;
    } else if (parameters.states[0].statetype == "GEqOE") {
        statetype = thames::constants::statetypes::GEQOE;
    } else if (parameters.states[0].statetype == "Keplerian") {
        statetype = thames::constants::statetypes::KEPLERIAN;
    } else {
        throw std::runtime_error("Unsupported state type provided");
    }

    // Propagator
    if (parameters.propagator.equations == "Cowell") {
        // Set up propagator
        auto propagator = thames::propagators::CowellPropagator<T>(mu, perturbation, factors);
        // Propagate
        states_propagated = propagator.propagate(tstart, tend, tstep, states, parameters.propagator, statetype);
    } else if (parameters.propagator.equations == "GEqOE") {
        // Set up propagator
        auto propagator = thames::propagators::GEqOEPropagator<T>(mu, perturbation, factors);
        // Propagate
        states_propagated = propagator.propagate(tstart, tend, tstep, states, parameters.propagator, statetype);        
    } else {
        throw std::runtime_error("Unsupported propagator requested");
    }

    // Declare output structure
    thames::settings::Parameters<T> parameters_output(parameters);

    // Set flag
    parameters_output.metadata.isInputFile = false;

    // Update output states
    thames::settings::StateParameters<T> state_output;
    state_output.datetime = tend;
    state_output.states = states_propagated;
    state_output.statetype = parameters.states[0].statetype;
    parameters_output.states.push_back(state_output);

    // Return parameters
    return parameters_output;
}

template<class T, template <class> class P>
thames::settings::Parameters<T> propagate(const thames::settings::Parameters<T>& parameters) {
    // Load constants
    T J2 = thames::constants::earth::J2;
    T mu = thames::constants::earth::mu;
    T radius = thames::constants::earth::radius;
    T w = thames::constants::earth::w;

    // Declare factors
    auto factors = std::make_shared<thames::conversions::dimensional::DimensionalFactors<T>>();

    // Set up perturbations
    auto perturbation = std::make_shared<thames::perturbations::perturbationcombiner::PerturbationCombinerPolynomial<T, P>>(factors);

    // Set up atmosphere model
    if (parameters.perturbation.atmosphere.isEnabled) {
        // Select atmosphere model
        if (parameters.perturbation.atmosphere.model == "USSA76") {
            auto atmospheremodel = thames::perturbations::atmosphere::models::USSA76;
            auto atmosphereperturbation = std::make_shared<thames::perturbations::atmosphere::drag::DragPolynomial<T, P>>(radius, w, parameters.spacecraft.Cd, parameters.spacecraft.dragArea, parameters.spacecraft.mass, atmospheremodel, factors);
            perturbation->add_model(atmosphereperturbation);
        } else if (parameters.perturbation.atmosphere.model == "Wertz") {
            auto atmospheremodel = thames::perturbations::atmosphere::models::WERTZ;
            auto atmosphereperturbation = std::make_shared<thames::perturbations::atmosphere::drag::DragPolynomial<T, P>>(radius, w, parameters.spacecraft.Cd, parameters.spacecraft.dragArea, parameters.spacecraft.mass, atmospheremodel, factors);
            perturbation->add_model(atmosphereperturbation);
        } else {
            throw std::runtime_error("Unsupported atmosphere model requested");
        }        
    }

    // Set up geopotential model
    if (parameters.perturbation.geopotential.isEnabled) {
        // Select geopotential model
        if (parameters.perturbation.geopotential.model == "J2") {
            auto geopotentialmodel = std::make_shared<thames::perturbations::geopotential::J2Polynomial<T, P>>(mu, J2, radius, factors);
            perturbation->add_model(geopotentialmodel);
        } else {
            throw std::runtime_error("Unsupported geopotential model requested");
        }
    }

    // Import states
    std::vector<std::vector<T>> states = parameters.states[0].states;
    std::vector<std::vector<T>> states_propagated;
    T tstart = parameters.propagator.startTime;
    T tend = parameters.propagator.endTime;
    T tstep = parameters.propagator.timeStep;

    // Import polynomial parameters
    unsigned int degree = parameters.polynomial.maxDegree;

    // Import state type
    thames::constants::statetypes::StateTypes statetype;
    if (parameters.states[0].statetype == "Cartesian") {
        statetype = thames::constants::statetypes::CARTESIAN;
    } else if (parameters.states[0].statetype == "GEqOE") {
        statetype = thames::constants::statetypes::GEQOE;
    } else if (parameters.states[0].statetype == "Keplerian") {
        statetype = thames::constants::statetypes::KEPLERIAN;
    } else {
        throw std::runtime_error("Unsupported state type provided");
    }

    // Propagator
    if (parameters.propagator.equations == "Cowell") {
        // Set up propagator
        auto propagator = thames::propagators::CowellPropagatorPolynomial<T, P>(mu, perturbation, factors);
        // Propagate
        states_propagated = propagator.propagate(tstart, tend, tstep, states, parameters.propagator, statetype, degree);
    } else if (parameters.propagator.equations == "GEqOE") {
        // Set up propagator
        auto propagator = thames::propagators::GEqOEPropagatorPolynomial<T, P>(mu, perturbation, factors);
        // Propagate
        states_propagated = propagator.propagate(tstart, tend, tstep, states, parameters.propagator, statetype, degree);        
    } else {
        throw std::runtime_error("Unsupported propagator requested");
    }

    // Declare output structure
    thames::settings::Parameters<T> parameters_output(parameters);

    // Set flag
    parameters_output.metadata.isInputFile = false;

    // Update output states
    thames::settings::StateParameters<T> state_output;
    state_output.datetime = tend;
    state_output.states = states_propagated;
    state_output.statetype = parameters.states[0].statetype;
    parameters_output.states.push_back(state_output);

    // Return parameters
    return parameters_output;
}

int main(int argc, char **argv) {
    // Throw error if incorrect number of arguments are provided
    if (argc != 3)
        throw std::runtime_error("Incorrect number of arguments provided");

    // Store filepaths as strings
    std::string filepathin(argv[1]), filepathout(argv[2]);

    // Load input file
    thames::settings::Parameters<double> parameters, parameters_output;
    thames::io::json::load(filepathin, parameters);

    // Throw error if polynomial propagation requested with version of THAMES not compiled with SMART-UQ support
    #ifndef THAMES_USE_SMARTUQ
    if (parameters.polynomial.isEnabled)
        throw std::runtime_error("Polynomial propagation requested for version of THAMES not compiled with SMART-UQ");
    #endif

    // Check valid input parameters
    if (!parameters.metadata.isInputFile)
        throw std::runtime_error("Input file flag set to false");
    if (parameters.states.size() != 1)
        throw std::runtime_error("Non-singular set of input states provided");
    if (parameters.states[0].datetime != parameters.propagator.startTime)
        throw std::runtime_error("Inconsistent start times provided");

    // Start timer for propagation
    auto start_propagation = std::chrono::high_resolution_clock::now();

    // Propagate
    if (parameters.polynomial.isEnabled) {
        if (parameters.polynomial.type == "Taylor") {
            parameters_output = propagate<double, smartuq::polynomial::taylor_polynomial>(parameters);
        } else if (parameters.polynomial.type == "Chebyshev") {
            parameters_output = propagate<double, smartuq::polynomial::chebyshev_polynomial>(parameters);
        } else {
            throw std::runtime_error("Unsupported polynomial type requested");
        }
    } else {
        parameters_output = propagate<double>(parameters);
    }

    // Start timer for propagation
    auto end_propagation = std::chrono::high_resolution_clock::now();

    // Store propagation statistics
    std::chrono::duration<double> elapsed_propagation = end_propagation - start_propagation;
    parameters_output.statistics.propagationTime = elapsed_propagation.count();

    // Output file
    thames::io::json::save(filepathout, parameters_output);

    return 0;
}