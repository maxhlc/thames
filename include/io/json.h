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

#ifndef THAMES_IO_JSON
#define THAMES_IO_JSON

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

namespace thames::io::json {

    /**
     * @brief Structure to store JSON file metadata
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-24
     * 
     */
    struct Metadata {
        /// File name
        std::string name;

        /// File description
        std::string description;

        /// File creation date
        std::string datetimeCreated;

        /// File modification date
        std::string datetimeModified;

        /// Input flag
        bool isInputFile;

        // Macro to generate boilerplate to/from JSON
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(Metadata, name, description, datetimeCreated, datetimeModified, isInputFile)
    };

    /**
     * @brief Structure to store spacecraft parameters
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-24
     * 
     * @tparam T Numeric type
     */
    template<class T>
    struct SpacecraftParameters {
        /// Spacecraft mass [kg]
        T mass;

        /// Spacecraft drag area [km^2]
        T dragArea;

        // Macro to generate boilerplate to/from JSON
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(SpacecraftParameters, mass, dragArea)
    };

    /**
     * @brief Structure to store geopotential model parameters
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-24
     * 
     */
    struct GeopotentialPerturbationParameters {
        /// Enabled flag
        bool isEnabled;

        /// Geopotential model
        std::string model;

        /// Maximum model order
        unsigned int maxOrder;

        /// Maximum model degree
        unsigned int maxDegree;

        // Macro to generate boilerplate to/from JSON
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(GeopotentialPerturbationParameters, isEnabled, model, maxOrder, maxDegree)
    };

    /**
     * @brief Structure to store atmospheric model parameters
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-24
     * 
     */
    struct AtmospherePerturbationParameters {
        /// Enabled flag
        bool isEnabled;

        /// Atmosphere model
        std::string model;

        // Macro to generate boilerplate to/from JSON
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(AtmospherePerturbationParameters, isEnabled, model)
    };

    /**
     * @brief Structure to store perturbation model parameters
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-24
     * 
     */
    struct PerturbationParameters {
        /// Geopotential model
        GeopotentialPerturbationParameters geopotential;

        /// Atmosphere model
        AtmospherePerturbationParameters atmosphere;

        // Macro to generate boilerplate to/from JSON
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(PerturbationParameters, geopotential, atmosphere)
    };

    /**
     * @brief Structure to store propagator parameters
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-26
     * 
     * @tparam T Numeric type
     */
    template<class T>
    struct PropagatorParameters {
        /// Propagation start time
        T startTime;

        /// Propagation end time
        T endTime;

        /// Propagation equations
        std::string equations;

        // Macro to generate boilerplate to/from JSON
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(PropagatorParameters, startTime, endTime, equations)
    };

    /**
     * @brief Structure to store integrator parameters
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-24
     * 
     * @tparam T Numeric type
     */
    template<class T>
    struct IntegratorParameters {
        /// Fixed- or variable-step flag
        bool isFixedStep;

        /// Intermediate output flag
        bool intermediateOutput;

        /// Fixed-step timestep
        T timeStep;

        /// Variable-step absolute tolerance
        T absoluteTolerance;

        /// Variable-step relative tolerance
        T relativeTolerance;

        // Macro to generate boilerplate to/from JSON
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(IntegratorParameters, isFixedStep, intermediateOutput, timeStep, absoluteTolerance, relativeTolerance)
    };

    /**
     * @brief Structure to store polynomial parameters
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-24
     * 
     */
    struct PolynomialParameters {
        /// Enabled flag
        bool isEnabled;

        /// Polynomial type
        std::string type;

        /// Maximum polynomial degree
        unsigned int maxDegree;

        // Macro to generate boilerplate to/from JSON
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(PolynomialParameters, isEnabled, type, maxDegree)
    };

    /**
     * @brief Structure to store state parameters
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-26
     * 
     * @tparam T Numeric type
     */
    template<class T>
    struct StateParameters {
        /// State time
        T datetime;

        /// State vectors
        std::vector<std::vector<T>> states;

        /// State type
        std::string statetype;

        // Macro to generate boilerplate to/from JSON
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(StateParameters, datetime, states, statetype)
    };

    /**
     * @brief Structure to store all parameters
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-26
     * 
     * @tparam T Numeric type
     */
    template<class T>
    struct Parameters {
        /// File metadata
        Metadata metadata;

        /// Spacecraft parameters
        SpacecraftParameters<T> spacecraft;

        /// Perturbation parameters
        PerturbationParameters perturbation;

        /// Propagator parameters
        PropagatorParameters<T> propagator;

        /// Integrator parameters
        IntegratorParameters<T> integrator;

        /// Polynomial parameters
        PolynomialParameters polynomial;

        /// State parameters
        std::vector<StateParameters<T>> states;

        // Macro to generate boilerplate to/from JSON
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(Parameters, metadata, spacecraft, perturbation, propagator, integrator, polynomial, states)
    };

    /**
     * @brief Load parameters from JSON
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-26
     * 
     * @tparam T Numeric type
     * @param[in] filepath Input file path
     * @param[out] parameters Parameters
     */
    template<class T>
    void load(const std::string& filepath, Parameters<T>& parameters);

    /**
     * @brief Save parameters to JSON
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-26
     * 
     * @tparam T Numeric type
     * @param[in] filepath Input file path
     * @param[in] parameters Parameters
     */
    template<class T>
    void save(const std::string& filepath, const Parameters<T>& parameters);

}

#endif