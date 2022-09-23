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

#ifndef THAMES_PERTURBATIONS_ATMOSPHERE_WERTZP1
#define THAMES_PERTURBATIONS_ATMOSPHERE_WERTZP1

#include <vector>

#include "baseatmospheremodel.h"

namespace thames::perturbations::atmosphere::models {

    /**
     * @brief Approximate Wertz Exponential Atmospheric Model
     * 
     * Polynomial approximation using values from Vallado, 2013
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-08-24
     * 
     * @tparam Numeric type
     */
    template<class T>
    class WertzP1AtmosphereModel : public BaseAtmosphereModel<T> {

        private:

            // Polynomial approximation domain [km]
            const std::vector<T> m_domain = {250.0, 1000.0};

            // Polynomial coefficients
            const std::vector<T> m_coeff = {-29.41808788, -5.10257202};

        public:

            /**
             * @brief Construct a new Approximate Wertz Atmosphere Model object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-24
             */
            WertzP1AtmosphereModel();

            /**
             * @brief Destroy the Approximate Wertz Atmosphere Model object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-24
             */
            ~WertzP1AtmosphereModel();

            /**
             * @brief Calculate density using a polynomial approximation of the density profile
             * 
             * @param[in] alt Altitude [km]
             * @return T Atmospheric density [kg/m^3]
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-24
             */
            T density(T alt) const override;

    };

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Approximate Wertz Exponential Atmospheric Model (Polynomial)
     * 
     * Polynomial approximation using values from Vallado, 2013
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-08-24
     * 
     * @tparam T Numeric type
     * @tparam P Polynomial type
     */
    template<class T, template <class> class P>
    class WertzP1AtmosphereModelPolynomial : public BaseAtmosphereModelPolynomial<T, P> {
        
        private:

            // Polynomial approximation domain [km]
            const std::vector<T> m_domain = {250.0, 1000.0};

            // Polynomial coefficients
            const std::vector<T> m_coeff = {-29.41808788, -5.10257202};

        public:

            /**
             * @brief Construct a new Approximate Wertz Atmosphere Model Polynomial object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-24
             */
            WertzP1AtmosphereModelPolynomial();

            /**
             * @brief Destroy the Approximate Wertz Atmosphere Model Polynomial object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-24
             */
            ~WertzP1AtmosphereModelPolynomial();

            /**
             * @brief Calculate density using a polynomial approximation of the density profile
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-24
             * 
             * @param[in] alt Altitude [km]
             * @return P<T> Atmospheric density [kg/m^3]
             */
            P<T> density(P<T> alt) const override;
    };

    #endif

}

#endif