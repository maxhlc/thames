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
#include <string>

#include <nlohmann/json.hpp>

#include "../../include/io/json.h"

namespace thames::io::json {

    template<class T>
    void load(const std::string& filepath, Metadata& metadata, SpacecraftParameters<T>& spacecraft, PerturbationParameters& perturbation, PropagatorParameters& propagator, IntegratorParameters<T>& integrator, StateParameters<T>& state) {
        /// @todo Input checking

        // Open file stream
        std::ifstream filestream(filepath);

        // Load JSON
        nlohmann::json j;
        filestream >> j;

        // Load metadata
        metadata = j["metadata"].get<Metadata>();

        // Throw error if not an input file
        if (!metadata.isInputFile)
            throw std::runtime_error("Input file flag must be set to true");

        // Throw error if more than one set of points is provided in the file
        if (j["states"].size() > 1)
            throw std::runtime_error("Input file must have only one set of points");

        // Load spacecraft parameters
        spacecraft = j["spacecraft"].get<SpacecraftParameters<T>>();

        // Load perturbation parameters
        perturbation = j["perturbation"].get<PerturbationParameters>();

        // Load propagator parameters
        propagator = j["propagator"].get<PropagatorParameters>();

        // Load integrator parameters
        integrator = j["integrator"].get<IntegratorParameters<T>>();

        // Load state parameters
        state = j["states"][0].get<StateParameters<T>>();

        // Check consistency between start times
        if (propagator.startTime != state.datetime)
            throw std::runtime_error("Start times do not match");
    }
    template void load(const std::string&, Metadata&, SpacecraftParameters<double>&, PerturbationParameters&, PropagatorParameters&, IntegratorParameters<double>&, StateParameters<double>&);

    template<class T>
    void save(const std::string& filepath, const Metadata& metadata, const SpacecraftParameters<T>& spacecraft, const PerturbationParameters& perturbation, const PropagatorParameters& propagator, const IntegratorParameters<T>& integrator, const std::vector<StateParameters<T>>& states) {
        // Open file stream
        std::ofstream filestream(filepath);

        // Construct JSON object
        nlohmann::json j = {
            {"metadata", metadata},
            {"spacecraft", spacecraft},
            {"perturbation", perturbation},
            {"propagator", propagator},
            {"integrator", integrator},
            {"states", states}
        };

        // Output JSON object
        filestream << std::setw(4) << j;
    }
    template void save(const std::string&, const Metadata&, const SpacecraftParameters<double>&, const PerturbationParameters&, const PropagatorParameters& , const IntegratorParameters<double>&, const std::vector<StateParameters<double>>&); 

    #ifdef THAMES_USE_SMARTUQ

    template<class T>
    void load(const std::string& filepath, Metadata& metadata, SpacecraftParameters<T>& spacecraft, PerturbationParameters& perturbation, PropagatorParameters& propagator, IntegratorParameters<T>& integrator, PolynomialParameters& polynomial, StateParameters<T>& state) {
        /// @todo Input checking

        // Open file stream
        std::ifstream filestream(filepath);

        // Load JSON
        nlohmann::json j;
        filestream >> j;

        // Load metadata
        metadata = j["metadata"].get<Metadata>();

        // Throw error if not an input file
        if (!metadata.isInputFile)
            throw std::runtime_error("Input file flag must be set to true");

        // Throw error if more than one set of points is provided in the file
        if (j["states"].size() > 1)
            throw std::runtime_error("Input file must have only one set of points");

        // Load spacecraft parameters
        spacecraft = j["spacecraft"].get<SpacecraftParameters<T>>();

        // Load perturbation parameters
        perturbation = j["perturbation"].get<PerturbationParameters>();

        // Load propagator parameters
        propagator = j["propagator"].get<PropagatorParameters>();

        // Load integrator parameters
        integrator = j["integrator"].get<IntegratorParameters<T>>();

        // Load polynomial parameters
        polynomial = j["polynomial"].get<PolynomialParameters>();

        // Load state parameters
        state = j["states"][0].get<StateParameters<T>>();

        // Check consistency between start times
        if (propagator.startTime != state.datetime)
            throw std::runtime_error("Start times do not match");
    }
    template void load(const std::string&, Metadata&, SpacecraftParameters<double>&, PerturbationParameters&, PropagatorParameters&, IntegratorParameters<double>&, PolynomialParameters&, StateParameters<double>&);

    template<class T>
    void save(const std::string& filepath, const Metadata& metadata, const SpacecraftParameters<T>& spacecraft, const PerturbationParameters& perturbation, const PropagatorParameters& propagator, const IntegratorParameters<T>& integrator, const PolynomialParameters& polynomial, const std::vector<StateParameters<T>>& states) {
        // Open file stream
        std::ofstream filestream(filepath);

        // Construct JSON object
        nlohmann::json j = {
            {"metadata", metadata},
            {"spacecraft", spacecraft},
            {"perturbation", perturbation},
            {"propagator", propagator},
            {"integrator", integrator},
            {"polynomial", polynomial},
            {"states", states}
        };

        // Output JSON object
        filestream << std::setw(4) << j;
    }
    template void save(const std::string&, const Metadata&, const SpacecraftParameters<double>&, const PerturbationParameters&, const PropagatorParameters& , const IntegratorParameters<double>&, const PolynomialParameters&, const std::vector<StateParameters<double>>&);

    #endif 

}