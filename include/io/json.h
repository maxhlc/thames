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
     * @brief Struct to store JSON file metadata
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
     * @brief Struct to store spacecraft parameters
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
     * @brief Struct to store geopotential model parameters
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
     * @brief Struct to store atmospheric model parameters
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
     * @brief Struct to store perturbation model parameters
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
     * @brief Struct to store propagator parameters
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-24
     * 
     */
    struct PropagatorParameters {
        /// Propagation start time
        std::string startTime;

        /// Propagation end time
        std::string endTime;

        /// Propagation equations
        std::string equations;

        // Macro to generate boilerplate to/from JSON
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(PropagatorParameters, startTime, endTime, equations)
    };

    /**
     * @brief Struct to store integrator parameters
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

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Struct to store polynomial parameters
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

    #endif

    /**
     * @brief Struct to store state parameters
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-24
     * 
     * @tparam T Numeric type
     */
    template<class T>
    struct StateParameters {
        /// State time
        std::string datetime;

        /// State vectors
        std::vector<std::vector<T>> states;

        /// State type
        std::string statetype;

        // Macro to generate boilerplate to/from JSON
        NLOHMANN_DEFINE_TYPE_INTRUSIVE(StateParameters, datetime, states, statetype)
    };

    /**
     * @brief Load parameters from JSON
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-24
     * 
     * @tparam T Numeric type
     * @param[in] filepath Input file path
     * @param[out] metadata File metadata
     * @param[out] spacecraft Spacecraft parameters
     * @param[out] perturbation Perturbation parameters
     * @param[out] propagator Propagator parameters
     * @param[out] integrator Integrator parameters
     * @param[out] state State parameters
     */
    template<class T>
    void load(const std::string& filepath, Metadata& metadata, SpacecraftParameters<T>& spacecraft, PerturbationParameters& perturbation, PropagatorParameters& propagator, IntegratorParameters<T>& integrator, StateParameters<T>& state);

    /**
     * @brief Save parameters to JSON
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-24
     * 
     * @tparam T Numeric type
     * @param[in] filepath Input file path
     * @param[in] metadata File metadata
     * @param[in] spacecraft Spacecraft parameters
     * @param[in] perturbation Perturbation parameters
     * @param[in] propagator Propagator parameters
     * @param[in] integrator Integrator parameters
     * @param[in] states State parameters
     */
    template<class T>
    void save(const std::string& filepath, const Metadata& metadata, const SpacecraftParameters<T>& spacecraft, const PerturbationParameters& perturbation, const PropagatorParameters& propagator, const IntegratorParameters<T>& integrator, const std::vector<StateParameters<T>>& states);

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Load parameters from JSON
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-24
     * 
     * @tparam T Numeric type
     * @param[in] filepath Input file path
     * @param[out] metadata File metadata
     * @param[out] spacecraft Spacecraft parameters
     * @param[out] perturbation Perturbation parameters
     * @param[out] propagator Propagator parameters
     * @param[out] integrator Integrator parameters
     * @param[out] polynomial Polynomial parameters
     * @param[out] state State parameters
     */
    template<class T>
    void load(const std::string& filepath, Metadata& metadata, SpacecraftParameters<T>& spacecraft, PerturbationParameters& perturbation, PropagatorParameters& propagator, IntegratorParameters<T>& integrator, PolynomialParameters& polynomial, StateParameters<T>& state);

    /**
     * @brief Save parameters to JSON
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-24
     * 
     * @tparam T Numeric type
     * @param[in] filepath Input file path
     * @param[in] metadata File metadata
     * @param[in] spacecraft Spacecraft parameters
     * @param[in] perturbation Perturbation parameters
     * @param[in] propagator Propagator parameters
     * @param[in] integrator Integrator parameters
     * @param[in] polynomial Polynomial parameters
     * @param[in] states State parameters
     */
    template<class T>
    void save(const std::string& filepath, const Metadata& metadata, const SpacecraftParameters<T>& spacecraft, const PerturbationParameters& perturbation, const PropagatorParameters& propagator, const IntegratorParameters<T>& integrator, const PolynomialParameters& polynomial, const std::vector<StateParameters<T>>& states);

    #endif

}

#endif