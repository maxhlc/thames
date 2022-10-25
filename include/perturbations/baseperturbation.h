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

#ifndef THAMES_PERTURBATIONS_BASEPERTURBATION
#define THAMES_PERTURBATIONS_BASEPERTURBATION

#include <array>
#include <memory>
#include <vector>

#include "../conversions/dimensional.h"

namespace thames::perturbations::baseperturbation{

    using thames::conversions::dimensional::DimensionalFactors;
    
    ///////////
    // Reals //
    ///////////

    /**
     * @brief Class for the base perturbation.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-27
     * 
     * @tparam T Numeric type.
     */
    template<class T>
    class BasePerturbation{

        protected:

            /// Dimensional factors
            const std::shared_ptr<const DimensionalFactors<T>> m_factors;

            /// Flag for whether to expect non-dimensional inputs/outputs
            bool m_isNonDimensional = false;

        public:

            /**
             * @brief Construct a new Base Perturbation object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-27
             * 
             * @param[in] factors Factors for non-dimensionalisation.
             * 
             */
            BasePerturbation(const std::shared_ptr<const DimensionalFactors<T>> factors);

            /**
             * @brief Destroy the Base Perturbation object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-10
             * 
             */
            ~BasePerturbation();

            /**
             * @brief Get the non-dimensional flag.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-19
             * 
             * @return bool
             */
            bool get_nondimensional();

            /**
             * @brief Set the non-dimensional flag.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-20
             * 
             * @param[in] isNonDimensional Non-dimensional flag.
             */
            virtual void set_nondimensional(const bool isNonDimensional);

            /**
             * @brief Default total perturbing acceleration.
             * 
             * Returns zero total perturbing acceleration.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-25
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::vector<T> Total perturbing acceleration.
             */
            virtual std::vector<T> acceleration_total(const T& t, const std::vector<T>& R, const std::vector<T>& V) const;

            /**
             * @brief Default non-potential perturbing acceleration.
             * 
             * Returns zero non-potential perturbing acceleration.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-25
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::vector<T> Non-potential perturbing acceleration.
             */
            virtual std::vector<T> acceleration_nonpotential(const T& t, const std::vector<T>& R, const std::vector<T>& V) const;

            /**
             * @brief Default perturbing potential.
             * 
             * Returns zero potential.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-25
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return T Perturbing potential.
             */
            virtual T potential(const T& t, const std::vector<T>& R) const;

            /**
             * @brief Default time derivative of the perturbing potential.
             * 
             * Returns zero time derivative of the perturbing potential.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-25
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return T Time derivative of the perturbing potential.
             */
            virtual T potential_derivative(const T& t, const std::vector<T>& R, const std::vector<T>& V) const;

    };

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Class for the base perturbation for polynomial distributions.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-27
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     */
    template<class T, template<class> class P>
    class BasePerturbationPolynomial{

        protected:

            /// Dimensional factors
            const std::shared_ptr<const DimensionalFactors<T>> m_factors;

            /// Flag for whether to expect non-dimensional inputs/outputs
            bool m_isNonDimensional = false;

        public:

            /**
             * @brief Construct a new Base Perturbation Polynomial object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-27
             * 
             * @param[in] factors Factors for non-dimensionalisation.
             * 
             */
            BasePerturbationPolynomial(const std::shared_ptr<const DimensionalFactors<T>> factors);

            /**
             * @brief Destroy the Base Perturbation Polynomial object.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-27
             * 
             */
            ~BasePerturbationPolynomial();

            /**
             * @brief Get the non-dimensional flag.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-11
             * 
             * @return bool
             */
            virtual bool get_nondimensional();

            /**
             * @brief Set the non-dimensional flag.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-05-11
             * 
             * @param isNonDimensional Non-dimensional flag.
             */
            virtual void set_nondimensional(bool isNonDimensional);

            /**
             * @brief Default total perturbing acceleration.
             * 
             * Returns zero total perturbing acceleration.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-27
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::vector<P<T>> Total perturbing acceleration.
             */
            virtual std::vector<P<T>> acceleration_total(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const;

            /**
             * @brief Default non-potential perturbing acceleration.
             * 
             * Returns zero non-potential perturbing acceleration.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-27
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::vector<P<T>> Non-potential perturbing acceleration.
             */
            virtual std::vector<P<T>> acceleration_nonpotential(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const;

            /**
             * @brief Default perturbing potential.
             * 
             * Returns zero potential.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-27
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return P<T> Perturbing potential.
             */
            virtual P<T> potential(const T& t, const std::vector<P<T>>& R) const;

            /**
             * @brief Default time derivative of the perturbing potential.
             * 
             * Returns zero time derivative of the perturbing potential.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-27
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return P<T> Time derivative of the perturbing potential.
             */
            virtual P<T> potential_derivative(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const;

    };

    #endif

}

#endif