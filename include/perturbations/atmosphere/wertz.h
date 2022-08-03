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

#ifndef THAMES_PERTURBATIONS_ATMOSPHERE_WERTZ
#define THAMES_PERTURBATIONS_ATMOSPHERE_WERTZ

#include <vector>

#include "baseatmospheremodel.h"

namespace thames::perturbations::atmosphere::models {

    /**
     * @brief Wertz Exponential Atmospheric Model
     * 
     * Values from Vallado, 2013
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-08-02
     * 
     * @tparam Numeric type
     */
    template<class T>
    class WertzAtmosphereModel : public BaseAtmosphereModel<T> {

        private:

            /// Geometric altitudes [km]
            const std::vector<T> m_geo = {
                0,    25,  30,  40,  50,  60,   70,
                80,   90, 100, 110, 120, 130,  140,
                150, 180, 200, 250, 300, 350,  400,
                450, 500, 600, 700, 800, 900, 1000
            };

            /// Atmospheric densities [kg/m^3]
            const std::vector<T> m_rho = {
                1.225E+00, 3.899E-02, 1.774E-02, 3.972E-03, 1.057E-03, 3.206E-04, 8.770E-05,
                1.905E-05, 3.396E-06, 5.297E-07, 9.661E-08, 2.438E-08, 8.484E-09, 3.845E-09,
                2.070E-09, 5.464E-10, 2.789E-10, 7.248E-11, 2.418E-11, 9.518E-12, 3.725E-12,
                1.585E-12, 6.967E-13, 1.454E-13, 3.614E-14, 1.170E-14, 5.245E-15, 3.019E-15
            };

            /// Scale heights [km]
            const std::vector<T> m_scale = {
                7.249,   6.349,  6.682,  7.554,   8.382,   7.714,   6.549,
                5.799,   5.382,  5.877,  7.263,   9.473,  12.636,  16.149,
                22.523, 29.740, 37.105, 45.546,  53.628,  53.298,  58.515,
                60.828, 63.822, 71.835, 88.667, 124.640, 181.050, 268.000
            };

        public:

            /**
             * @brief Construct a new Wertz Atmosphere Model object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             */
            WertzAtmosphereModel();

            /**
             * @brief Destroy the Wertz Atmosphere Model object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             */
            ~WertzAtmosphereModel();

            /**
             * @brief Calculate density using exponetial interpolation
             * 
             * @param[in] alt Altitude [km]
             * @return T Atmospheric density [kg/m^3]
             */
            T density(T alt) const override;

    };

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Wertz Exponential Atmospheric Model (Polynomial)
     * 
     * Values from Vallado, 2013
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-08-02
     * 
     * @tparam T Numeric type
     * @tparam P Polynomial type
     */
    template<class T, template <class> class P>
    class WertzAtmosphereModelPolynomial : public BaseAtmosphereModelPolynomial<T, P> {
        
        private:

        public:

            /**
             * @brief Construct a new Wertz Atmosphere Model Polynomial object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             */
            WertzAtmosphereModelPolynomial();

            /**
             * @brief Destroy the Wertz Atmosphere Model Polynomial object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             */
            ~WertzAtmosphereModelPolynomial();

            /**
             * @brief Calculate density using a pre-computed polynomial representation of the density profile
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             * 
             * @param[in] alt Altitude [km]
             * @return P<T> Atmospheric density [kg/m^3]
             */
            P<T> density(P<T> alt) const override;
    };

    #endif

}

#endif