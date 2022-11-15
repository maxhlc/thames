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

#ifndef THAMES_PERTURBATIONS_ATMOSPHERE_DRAG
#define THAMES_PERTURBATIONS_ATMOSPHERE_DRAG

#include "baseatmospheremodel.h"
#include "../baseperturbation.h"
#include "../../conversions/dimensional.h"

namespace thames::perturbations::atmosphere::drag {

    using thames::perturbations::atmosphere::models::BaseAtmosphereModel;
    using thames::perturbations::baseperturbation::BasePerturbation;
    using thames::conversions::dimensional::DimensionalFactors;

    /**
     * @brief Class for the perturbation resulting from atmospheric drag.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-08-02
     * 
     * @tparam T Numeric type.
     */
    template<class T>
    class Drag : public BasePerturbation<T> {

        private:

            /// Dimensional factors
            using BasePerturbation<T>::m_factors;

            /// Non-dimensional flag
            using BasePerturbation<T>::m_isNonDimensional;

            /// Central body radius
            const T m_radius;

            /// Central body rotation rate
            const T m_w;

            /// Spacecraft drag coefficient
            const T m_Cd;

            /// Spacecraft drag area
            const T m_A;

            /// Spacecraft mass
            const T m_m;

            /// Atmosphere model
            const std::shared_ptr<const BaseAtmosphereModel<T>> m_model;

        public:

            /**
             * @brief Construct a new Drag object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             * 
             * @param[in] radius Central body radius.
             * @param[in] w Central body rotation rate.
             * @param[in] Cd Spacecraft drag coefficient.
             * @param[in] A Spacecraft drag area.
             * @param[in] m Spacecraft mass.
             * @param[in] model Atmosphere model.
             * @param[in] factors Dimensional factors.
             */
            Drag(const T& radius, const T& w, const T& Cd, const T& A, const T& m, const std::shared_ptr<const BaseAtmosphereModel<T>> model, const std::shared_ptr<const DimensionalFactors<T>> factors);

            /**
             * @brief Destroy the Drag object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-25
             * 
             */
            ~Drag();

            /**
             * @brief Calculate perturbing acceleration resulting from drag. 
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-17
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::vector<T> Total perturbing acceleration due to drag.
             */
            std::vector<T> acceleration_total(const T& t, const std::vector<T>& R, const std::vector<T>& V) const override;

            /**
             * @brief Calculate perturbing acceleration resulting from drag. 
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::vector<T> Non-potential perturbing acceleration due to drag.
             */
            std::vector<T> acceleration_nonpotential(const T& t, const std::vector<T>& R, const std::vector<T>& V) const override;

    };

    #ifdef THAMES_USE_SMARTUQ

    using thames::perturbations::atmosphere::models::BaseAtmosphereModelPolynomial;
    using thames::perturbations::baseperturbation::BasePerturbationPolynomial;

    /**
     * @brief Class for the perturbation resulting from atmospheric drag.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-08-02
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     */
    template<class T, template <class> class P>
    class DragPolynomial : public BasePerturbationPolynomial<T, P> {

        private:

            /// Dimensional factors
            using BasePerturbationPolynomial<T, P>::m_factors;

            /// Non-dimensional flag
            using BasePerturbationPolynomial<T, P>::m_isNonDimensional;

            /// Central body radius
            const T m_radius;

            /// Central body rotation rate
            const T m_w;

            /// Spacecraft drag coefficient
            const T m_Cd;

            /// Spacecraft drag area
            const T m_A;

            /// Spacecraft mass
            const T m_m;

            /// Atmosphere model
            const std::shared_ptr<const BaseAtmosphereModelPolynomial<T, P>> m_model;

        public:

            /**
             * @brief Construct a new Drag object for use with polynomials.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             * 
             * @param[in] radius Central body radius.
             * @param[in] w Central body rotation rate.
             * @param[in] Cd Spacecraft drag coefficient.
             * @param[in] A Spacecraft drag area.
             * @param[in] m Spacecraft mass.
             * @param[in] model Atmosphere model.
             * @param[in] factors Dimensional factors.
             */
            DragPolynomial(const T& radius, const T& w, const T& Cd, const T& A, const T& m, const std::shared_ptr<const BaseAtmosphereModelPolynomial<T, P>> model, const std::shared_ptr<const DimensionalFactors<T>> factors);

            /**
             * @brief Destroy the Drag object for use with polynomials.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-17
             * 
             */
            ~DragPolynomial();

            /**
             * @brief Calculate perturbing acceleration resulting from drag. 
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-17
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::vector<P<T>> Total perturbing acceleration due to drag.
             */
            std::vector<P<T>> acceleration_total(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const override;

            /**
             * @brief Calculate perturbing acceleration resulting from drag. 
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-08-02
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::vector<P<T>> Non-potential perturbing acceleration due to drag.
             */
            std::vector<P<T>> acceleration_nonpotential(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const override;

    };

    #endif

}

#endif