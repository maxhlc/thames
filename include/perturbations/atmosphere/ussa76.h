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

#ifndef THAMES_PERTURBATIONS_ATMOSPHERE_USSA76
#define THAMES_PERTURBATIONS_ATMOSPHERE_USSA76

#include <vector>

#include "baseatmospheremodel.h"

namespace thames::perturbations::atmosphere::models {

    /**
     * @brief US Standard Atmosphere 1976
     * 
     * Values from Curtis, 2014
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-08-02
     * 
     * @tparam T Numeric type
     */
    template<class T>
    class USSA76AtmosphereModel : public BaseAtmosphereModel<T> {

        private:

            /// Geometric altitudes [km]
            const std::vector<T> m_geo = {
                  0,  25,  30,  40,  50,  60,   70,
                 80,  90, 100, 110, 120, 130,  140,
                150, 180, 200, 250, 300, 350,  400,
                450, 500, 600, 700, 800, 900, 1000
            };

            /// Atmospheric densities [kg/m^3]
            const std::vector<T> m_rho = {
                    1.225,  4.008e-2,  1.841e-2,  3.996e-3,  1.027e-3,  3.097e-4,  8.283e-5,
                 1.846e-5,  3.416e-6,  5.606e-7,  9.708e-8,  2.222e-8,  8.152e-9,  3.831e-9,
                 2.076e-9, 5.194e-10, 2.541e-10, 6.073e-11, 1.916e-11, 7.014e-12, 2.803e-12,
                1.184e-12, 5.215e-13, 1.137e-13, 3.070e-14, 1.136e-14, 5.759e-15, 3.561e-15
            };

            /// Scale heights [km]
            const std::vector<T> m_scale = {
                 7.310,  6.427,  6.546,   7.360,   8.342,   7.583,  6.661,
                 5.927,  5.533,  5.703,   6.782,   9.973,  13.243, 16.322,
                21.652, 27.974, 34.934,  43.342,  49.755,  54.513, 58.019,
                60.980, 65.654, 76.377, 100.587, 147.203, 208.020
            };

        public:

            /**
             * @brief Construct a new USSA76 Atmosphere Model object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             */
            USSA76AtmosphereModel();

            /**
             * @brief Destroy the USSA76 Atmosphere Model object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             */
            ~USSA76AtmosphereModel();

            /**
             * @brief Calculate density using exponential interpolation
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
     * @brief US Standard Atmosphere 1976 (Polynomial)
     * 
     * Values from Curtis, 2014
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-08-22
     * 
     * @tparam T Numeric type
     * @tparam P Polynomial type
     */
    template<class T, template <class> class P>
    class USSA76AtmosphereModelPolynomial : public BaseAtmosphereModelPolynomial<T, P> {

        private:

        public:

            /**
             * @brief Construct a new USSA76 Atmosphere Model Polynomial object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             */
            USSA76AtmosphereModelPolynomial();

            /**
             * @brief Destroy the USSA76 Atmosphere Model Polynomial object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             */
            ~USSA76AtmosphereModelPolynomial();

            /**
             * @brief Calculate density using exponential interpolation
             * 
             * @warning Not implemented
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-22
             * 
             * @param[in] alt Altitude [km]
             * @return P<T> Atmospheric density [kg/m^3]
             */
            P<T> density(P<T> alt) const override;
    };

    #endif

}

#endif