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

#ifndef THAMES_PERTURBATIONS_ATMOSPHERE_WERTZE
#define THAMES_PERTURBATIONS_ATMOSPHERE_WERTZE

#include <vector>

#include "baseatmospheremodel.h"

namespace thames::perturbations::atmosphere::models {

    /**
     * @brief Approximate Wertz Exponential Atmospheric Model
     * 
     * Exponential approximation using values from Vallado, 2013
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-08-23
     * 
     * @tparam Numeric type
     */
    template<class T>
    class WertzEAtmosphereModel : public BaseAtmosphereModel<T> {

        private:

            /// Exponential coefficients
            const T m_a = 1.0;
            const T m_b = -0.02068349;
            const T m_c = -18.19556173;
            const T m_d = 0.0;

        public:

            /**
             * @brief Construct a new Approximate Wertz Atmosphere Model object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-23
             */
            WertzEAtmosphereModel();

            /**
             * @brief Destroy the Approximate Wertz Atmosphere Model object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-23
             */
            ~WertzEAtmosphereModel();

            /**
             * @brief Calculate density using an exponential approximation of the density profile
             * 
             * @param[in] alt Altitude [km]
             * @return T Atmospheric density [kg/m^3]
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-23
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
     * Exponential approximation using values from Vallado, 2013
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-08-23
     * 
     * @tparam T Numeric type
     * @tparam P Polynomial type
     */
    template<class T, template <class> class P>
    class WertzEAtmosphereModelPolynomial : public BaseAtmosphereModelPolynomial<T, P> {
        
        private:

            /// Exponential coefficients
            const T m_a = 1.0;
            const T m_b = -0.02068349;
            const T m_c = -18.19556173;
            const T m_d = 0.0;

        public:

            /**
             * @brief Construct a new Approximate Wertz Atmosphere Model Polynomial object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-23
             */
            WertzEAtmosphereModelPolynomial();

            /**
             * @brief Destroy the Approximate Wertz Atmosphere Model Polynomial object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-23
             */
            ~WertzEAtmosphereModelPolynomial();

            /**
             * @brief Calculate density using an exponential approximation of the density profile
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-23
             * 
             * @param[in] alt Altitude [km]
             * @return P<T> Atmospheric density [kg/m^3]
             */
            P<T> density(P<T> alt) const override;
    };

    #endif

}

#endif