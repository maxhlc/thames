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

#ifndef THAMES_PERTURBATIONS_ATMOSPHERE_MODELS
#define THAMES_PERTURBATIONS_ATMOSPHERE_MODELS

#include <vector>

#include "ussa76.h"

namespace thames::perturbations::atmosphere::models {

    /// Enumerate for atmosphere models
    enum AtmosphereModels {
        USSA76
    };

    ///////////
    // Reals //
    ///////////

    /**
     * @brief Retrieve geometric altitudes
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-17
     * 
     * @tparam T Numeric type
     * @param[in] model Atmosphere model
     * @return std::vector<T> Geometric altitudes
     */
    template<class T>
    std::vector<T> get_geo(const AtmosphereModels& model);

    /**
     * @brief Retrieve atmospheric densities
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-17
     * 
     * @tparam T Numeric type
     * @param[in] model Atmosphere model
     * @return std::vector<T> Atmospheric densities
     */
    template<class T>
    std::vector<T> get_rho(const AtmosphereModels& model);

    /**
     * @brief Retrieve scale heights
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-17
     * 
     * @tparam T Numeric type
     * @param[in] model Atmosphere model
     * @return std::vector<T> Scale heights
     */
    template<class T>
    std::vector<T> get_scale(const AtmosphereModels& model);

    /**
     * @brief Atmosphere Model
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-17
     * 
     * @tparam T Numeric type
     */
    template<class T>
    class AtmosphereModel {

        private:

            /// Geometric altitudes
            const std::vector<T> m_geo;

            /// Atmospheric densities
            const std::vector<T> m_rho;

            /// Scale heights
            const std::vector<T> m_scale;

        public:

            /**
             * @brief Construct a new Atmosphere Model object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-17
             * 
             * @param[in] geo Geometric altitudes
             * @param[in] rho Atmospheric densities
             * @param[in] scale Scale heights
             */
            AtmosphereModel(const std::vector<T>& geo, const std::vector<T>& rho, const std::vector<T>& scale);

            /**
             * @brief Construct a new Atmosphere Model object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-17
             * 
             * @param[in] model Atmosphere model
             */
            AtmosphereModel(const AtmosphereModels& model);

            /**
             * @brief Destroy the Atmosphere Model object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-17
             * 
             */
            ~AtmosphereModel();

            /**
             * @brief Calculate density
             * 
             * Modified from Curtis, 2014
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-17
             * 
             * @param[in] alt Altitude
             * @return T Numeric type
             */
            T density(T alt) const;

    };

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Atmosphere Model
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-17
     * 
     * @tparam T Numeric type
     */
    template<class T, template <class> class P>
    class AtmosphereModelPolynomial {

        private:

            /// Geometric altitudes
            const std::vector<T> m_geo;

            /// Atmospheric densities
            const std::vector<T> m_rho;

            /// Scale heights
            const std::vector<T> m_scale;

        public:

            /**
             * @brief Construct a new Atmosphere Model object for use with polynomials
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-17
             * 
             * @param[in] geo Geometric altitudes
             * @param[in] rho Atmospheric densities
             * @param[in] scale Scale heights
             */
            AtmosphereModelPolynomial(const std::vector<T>& geo, const std::vector<T>& rho, const std::vector<T>& scale);

            /**
             * @brief Construct a new Atmosphere Model object for use with polynomials
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-17
             * 
             * @param[in] model Atmosphere model
             */
            AtmosphereModelPolynomial(const AtmosphereModels& model);

            /**
             * @brief Destroy the Atmosphere Model object for use with polynomials
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-17
             * 
             */
            ~AtmosphereModelPolynomial();

            /**
             * @brief Calculate density
             * 
             * Modified from Curtis, 2014
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-17
             * 
             * @param[in] alt Altitude
             * @return T Numeric type
             */
            P<T> density(P<T> alt) const;

    };

    #endif

}

#endif