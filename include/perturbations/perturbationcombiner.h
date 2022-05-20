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
     * @date 2022-05-20
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
             * @date 2022-05-20
             * 
             * @param[in] factors Dimensional factors
             */
            PerturbationCombiner(const DimensionalFactors<T>* const factors);

            /**
             * @brief Construct a new Perturbation Combiner object
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] models Perturbation models
             * @param[in] factors Dimensional factors
             */
            PerturbationCombiner(const std::vector<std::shared_ptr<BasePerturbation<T>>>& models, const DimensionalFactors<T>* const factors);

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

}

#endif