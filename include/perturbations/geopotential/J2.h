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

#ifndef THAMES_PERTURBATIONS_GEOPOTENTIAL_J2
#define THAMES_PERTURBATIONS_GEOPOTENTIAL_J2

#include <array>
#include <vector>

#include "../baseperturbation.h"
#include "../../conversions/dimensional.h"

namespace thames::perturbations::geopotential{

    using thames::perturbations::baseperturbation::BasePerturbation;
    using thames::conversions::dimensional::DimensionalFactors;

    ///////////
    // Reals //
    ///////////

    /**
     * @brief Class for the perturbation resulting from the J2-term.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-11
     * 
     * @tparam T Numeric type.
     */
    template<class T>
    class J2 : public BasePerturbation<T> {

        private:

            /// Dimensional factors
            using BasePerturbation<T>::m_factors;

            /// Non-dimensional flag
            using BasePerturbation<T>::m_isNonDimensional;

            /// Central body gravitational parameter
            const T m_mu;

            /// Central body J2-term
            const T m_J2;

            /// Central body radius
            const T m_radius;

        public:

            /**
             * @brief Construct a new J2 object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-11
             * 
             * @param[in] mu Central body gravitational parameter.
             * @param[in] J2 Central body J2-term.
             * @param[in] radius Central body radius.
             * @param[in] factors Dimensional factors.
             */
            J2(const T& mu, const T& J2, const T& radius, const DimensionalFactors<T>* factors);

            /**
             * @brief Destroy the J2 object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-25
             * 
             */
            ~J2();

            ////////////
            // Arrays //
            ////////////

            /**
             * @brief Calculate perturbing acceleration resulting from the J2-term. 
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-25
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::array<T, 3> Total perturbing acceleration due to the J2-term.
             */
            std::array<T, 3> acceleration_total(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const override;

            /**
             * @brief Calculate perturbing potential resulting from the J2-term. 
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-25
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return T Perturbing potential due to the J2-term.
             */
            T potential(const T& t, const std::array<T, 3>& R) const override;

            /////////////
            // Vectors //
            /////////////

            /**
             * @brief Calculate perturbing acceleration resulting from the J2-term. 
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-25
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::vector<T> Total perturbing acceleration due to the J2-term.
             */
            std::vector<T> acceleration_total(const T& t, const std::vector<T>& R, const std::vector<T>& V) const override;

            /**
             * @brief Calculate perturbing potential resulting from the J2-term. 
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-25
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return T Perturbing potential due to the J2-term.
             */
            T potential(const T& t, const std::vector<T>& R) const override;

    };

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using thames::perturbations::baseperturbation::BasePerturbationPolynomial;

    /**
     * @brief Class for the perturbation resulting from the J2-term.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-11
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     */
    template<class T, template<class> class P>
    class J2Polynomial : public BasePerturbationPolynomial<T, P> {
        
        private:

            /// Dimensional factors
            using BasePerturbationPolynomial<T, P>::m_factors;

            /// Non-dimensional flag
            using BasePerturbationPolynomial<T, P>::m_isNonDimensional;

            /// Central body gravitational parameter
            const T m_mu;       

            /// Central body J2-term
            const T m_J2;

            /// Central body radius
            const T m_radius;

        public:

            /**
             * @brief Construct a new J2Polynomial object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-11
             * 
             * @param[in] mu Central body gravitational parameter.
             * @param[in] J2 Central body J2-term.
             * @param[in] radius Central body radius.
             * @param[in] factors Dimensional factors.
             */
            J2Polynomial(const T& mu, const T& J2, const T& radius, const DimensionalFactors<T>* factors);

            /**
             * @brief Destroy the J2Polynomial object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-27
             * 
             */
            ~J2Polynomial();

            /**
             * @brief Calculate perturbing acceleration resulting from the J2-term. 
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-27
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::vector<P<T>> Total perturbing acceleration due to the J2-term.
             */
            std::vector<P<T>> acceleration_total(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const override;

            /**
             * @brief Calculate perturbing potential resulting from the J2-term. 
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-27
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return P<T> Perturbing potential due to the J2-term.
             */
            P<T> potential(const T& t, const std::vector<P<T>>& R) const override;

    };

    #endif

}

#endif