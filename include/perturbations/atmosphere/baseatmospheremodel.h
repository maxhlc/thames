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

namespace thames::perturbations::atmosphere::models {

    ///////////
    // Reals //
    ///////////

    /**
     * @brief Base Atmosphere Model
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-08-02
     * 
     * @tparam T Numeric type
     */
    template<class T>
    class BaseAtmosphereModel {

        private:

        public:

            /**
             * @brief Construct a new Base Atmosphere Model object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             * 
             */
            BaseAtmosphereModel();

            /**
             * @brief Destroy the Base Atmosphere Model object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             * 
             */
            ~BaseAtmosphereModel();

            /**
             * @brief Calculate density
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             * 
             * @param[in] alt Altitude
             * @return T Numeric type
             */
            virtual T density(T alt) const;

    };

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Base Atmosphere Model
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-08-02
     * 
     * @tparam T Numeric type
     * @tparam P Polynomial type
     */
    template<class T, template <class> class P>
    class BaseAtmosphereModelPolynomial {

        private:

        public:

            /**
             * @brief Construct a new Base Atmosphere Model object for use with polynomials
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             */
            BaseAtmosphereModelPolynomial();

            /**
             * @brief Destroy the Base Atmosphere Model object for use with polynomials
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             */
            ~BaseAtmosphereModelPolynomial();

            /**
             * @brief Calculate density
             * 
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             * 
             * @param[in] alt Altitude
             * @return T Numeric type
             */
            virtual P<T> density(P<T> alt) const;

    };

    #endif

}

#endif