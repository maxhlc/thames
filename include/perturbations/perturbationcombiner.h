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

#ifndef THAMES_PERTURBATIONS_PERTURBATIONCOMBINER
#define THAMES_PERTURBATIONS_PERTURBATIONCOMBINER

#include <memory>
#include <vector>

#include "baseperturbation.h"

namespace thames::perturbations::perturbationcombiner {

    using thames::conversions::dimensional::DimensionalFactors;
    using thames::perturbations::baseperturbation::BasePerturbation;
    
    /**
     * @brief Class to combine multiple perturbations
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-27
     * 
     * @tparam T Numeric type
     */
    template<class T>
    class PerturbationCombiner : public BasePerturbation<T> {

        private:

            /// Dimensional factors
            using BasePerturbation<T>::m_factors;

            /// Non-dimensional flag
            using BasePerturbation<T>::m_isNonDimensional;

            /// Underlying perturbation models
            std::vector<std::shared_ptr<BasePerturbation<T>>> m_models;

        public:

            /**
             * @brief Construct a new Perturbation Combiner object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-27
             * 
             * @param[in] factors Dimensional factors
             */
            PerturbationCombiner(const std::shared_ptr<const DimensionalFactors<T>> factors);

            /**
             * @brief Construct a new Perturbation Combiner object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-27
             * 
             * @param[in] models Perturbation models
             * @param[in] factors Dimensional factors
             */
            PerturbationCombiner(const std::vector<std::shared_ptr<BasePerturbation<T>>>& models, const std::shared_ptr<const DimensionalFactors<T>> factors);

            /**
             * @brief Destroy the Perturbation Combiner object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             */
            ~PerturbationCombiner();

            /**
             * @brief Set the non-dimensional flag
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] isNonDimensional Non-dimensional flag
             */
            void set_nondimensional(const bool isNonDimensional) override;

            /**
             * @brief Add perturbation model to combiner
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] model 
             */
            void add_model(const std::shared_ptr<BasePerturbation<T>>& model);

            ////////////
            // Arrays //
            ////////////

            /**
             * @brief Total perturbing acceleration
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] t Current physical time
             * @param[in] R Position vector
             * @param[in] V Velocity vector
             * @return std::array<T, 3> Total perturbing acceleration
             */
            std::array<T, 3> acceleration_total(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const override;

            /**
             * @brief Non-potential perturbing acceleration
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] t Current physical time
             * @param[in] R Position vector
             * @param[in] V Velocity vector
             * @return std::array<T, 3> Non-potential perturbing acceleration
             */
            std::array<T, 3> acceleration_nonpotential(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const override;

            /**
             * @brief Perturbing potential
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] t Current physical time
             * @param[in] R Position vector
             * @return T Perturbing potential
             */
            T potential(const T& t, const std::array<T, 3>& R) const override;

            /**
             * @brief Time derivative of the perturbing potential
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] t Current physical time
             * @param[in] R Position vector
             * @param[in] V Velocity vector
             * @return T Time derivative of the perturbing potential
             */
            T potential_derivative(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const override;

            /////////////
            // Vectors //
            /////////////

            /**
             * @brief Total perturbing acceleration
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] t Current physical time
             * @param[in] R Position vector
             * @param[in] V Velocity vector
             * @return std::vector<T> Total perturbing acceleration
             */
            std::vector<T> acceleration_total(const T& t, const std::vector<T>& R, const std::vector<T>& V) const override;

            /**
             * @brief Non-potential perturbing acceleration
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] t Current physical time
             * @param[in] R Position vector
             * @param[in] V Velocity vector
             * @return std::vector<T> Non-potential perturbing acceleration
             */
            std::vector<T> acceleration_nonpotential(const T& t, const std::vector<T>& R, const std::vector<T>& V) const override;

            /**
             * @brief Perturbing potential
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] t Current physical time
             * @param[in] R Position vector
             * @return T Perturbing potential
             */
            T potential(const T& t, const std::vector<T>& R) const override;

            /**
             * @brief Time derivative of the perturbing potential
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] t Current physical time
             * @param[in] R Position vector
             * @param[in] V Velocity vector
             * @return T Time derivative of the perturbing potential
             */
            T potential_derivative(const T& t, const std::vector<T>& R, const std::vector<T>& V) const override;
        
    };

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using thames::perturbations::baseperturbation::BasePerturbationPolynomial;

    /**
     * @brief Class to combine multiple polynomial perturbations
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-27
     * 
     * @tparam T Numeric type
     * @tparam P Polynomial type
     */
    template<class T, template <class> class P>
    class PerturbationCombinerPolynomial : public BasePerturbationPolynomial<T, P> {

        private:

            /// Dimensional factors
            using BasePerturbationPolynomial<T, P>::m_factors;

            /// Non-dimensional flag
            using BasePerturbationPolynomial<T, P>::m_isNonDimensional;

            /// Underlying perturbation models
            std::vector<std::shared_ptr<BasePerturbationPolynomial<T, P>>> m_models;

        public:

            /**
             * @brief Construct a new Perturbation Combiner object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-27
             * 
             * @param[in] factors Dimensional factors
             */
            PerturbationCombinerPolynomial(const std::shared_ptr<const DimensionalFactors<T>> factors);

            /**
             * @brief Construct a new Perturbation Combiner object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-27
             * 
             * @param[in] models Perturbation models
             * @param[in] factors Dimensional factors
             */
            PerturbationCombinerPolynomial(const std::vector<std::shared_ptr<BasePerturbationPolynomial<T, P>>>& models, const std::shared_ptr<const DimensionalFactors<T>> factors);

            /**
             * @brief Destroy the Perturbation Combiner object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             */
            ~PerturbationCombinerPolynomial();

            /**
             * @brief Set the non-dimensional flag
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] isNonDimensional Non-dimensional flag
             */
            void set_nondimensional(const bool isNonDimensional) override;

            /**
             * @brief Add perturbation model to combiner
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] model 
             */
            void add_model(const std::shared_ptr<BasePerturbationPolynomial<T, P>>& model);

            /**
             * @brief Total perturbing acceleration
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] t Current physical time
             * @param[in] R Position vector
             * @param[in] V Velocity vector
             * @return std::vector<P<T>> Total perturbing acceleration
             */
            std::vector<P<T>> acceleration_total(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const override;

            /**
             * @brief Non-potential perturbing acceleration
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] t Current physical time
             * @param[in] R Position vector
             * @param[in] V Velocity vector
             * @return std::vector<P<T>> Non-potential perturbing acceleration
             */
            std::vector<P<T>> acceleration_nonpotential(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const override;

            /**
             * @brief Perturbing potential
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] t Current physical time
             * @param[in] R Position vector
             * @return P<T> Perturbing potential
             */
            P<T> potential(const T& t, const std::vector<P<T>>& R) const override;

            /**
             * @brief Time derivative of the perturbing potential
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] t Current physical time
             * @param[in] R Position vector
             * @param[in] V Velocity vector
             * @return P<T> Time derivative of the perturbing potential
             */
            P<T> potential_derivative(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const override;
        
    };

    #endif

}

#endif