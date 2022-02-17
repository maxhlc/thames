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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../src/thames.h"

int main(int argc, char **argv){
    // Set constants
    double mu = thames::constants::earth::mu;

    // Store filepaths as strings
    std::string filepathin(argv[1]), filepathout(argv[2]);

    // Total times
    double total = 0;

    // Load sample states
    auto start = std::chrono::high_resolution_clock::now();
    double tstart, tend;
    int scid;
    std::vector<std::vector<double>> states;
    thames::io::point::load(filepathin, tstart, tend, scid, states);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Input time: " << duration.count() << " us \n";
    total += duration.count();

    // Declare propagator and perturbations
    thames::perturbations::baseperturbation::BasePerturbation<double> perturbation;
    thames::propagators::CowellPropagator<double> propagator(mu, &perturbation);

    // Declare vector for propagated states
    std::vector<std::vector<double>> states_propagated(states.size(), std::vector<double>(6));

    // Declare temporary vectors
    std::vector<double> state(6), state_propagated(6);

    // Iterate through samples
    start = std::chrono::high_resolution_clock::now();
    for(std::size_t ii=0; ii<states.size(); ii++){
        // Extract state
        state = states[ii];

        // Propagate state
        state_propagated = propagator.propagate(tstart, tend, 30, state);

        // Store propagated state
        states_propagated[ii] = state_propagated;
    }
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Propagation time: " << duration.count() << " us \n";
    total += duration.count();

    // Save propagated states
    start = std::chrono::high_resolution_clock::now();
    thames::io::point::save(filepathout, tstart, tend, scid, states_propagated);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Output time: " << duration.count() << " us \n";
    total += duration.count();

    // Print total time
    std::cout << "Total time: " << total << " us\n";

    // Return zero
    return 0;
}