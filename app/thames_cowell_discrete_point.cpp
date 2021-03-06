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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../include/thames.h"

int main(int argc, char **argv){
    // Set constants
    double mu = thames::constants::earth::mu;
    double radius = thames::constants::earth::radius;
    double J2 = thames::constants::earth::J2;
    double atol, rtol;

    // Store filepaths as strings
    std::string filepathin(argv[1]), filepathout(argv[2]);

    // Load sample states
    double tstart, tend, tstep = 30;
    int scid, degree;
    thames::constants::statetypes::StateTypes statetype;
    std::vector<std::vector<double>> states;
    thames::io::point::load(filepathin, tstart, tend, scid, statetype, degree, atol, rtol, states);

    // Create propagator options
    thames::settings::PropagatorParameters<double> options;
    options.absoluteTolerance = atol;
    options.relativeTolerance = rtol;

    // Calculate non-dimensionalisation factors based on first point
    auto factors = std::make_shared<thames::conversions::dimensional::DimensionalFactors<double>>();

    // Declare propagator and perturbations
    auto perturbation = std::make_shared<thames::perturbations::geopotential::J2<double>>(mu, J2, radius, factors);
    thames::propagators::CowellPropagator<double> propagator(mu, perturbation, factors);

    // Declare vector for propagated states
    std::vector<std::vector<double>> states_propagated(states.size(), std::vector<double>(6));

    // Declare temporary vectors
    std::vector<double> state(6), state_propagated(6);

    // Iterate through samples
    states_propagated = propagator.propagate(tstart, tend, tstep, states, options, statetype);

    // Save propagated states
    thames::io::point::save(filepathout, tstart, tend, scid, statetype, degree, atol, rtol, states_propagated);

    // Return zero
    return 0;
}