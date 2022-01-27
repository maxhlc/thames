#ifndef THAMES_PERTURBATIONS_BASEPERTURBATION
#define THAMES_PERTURBATIONS_BASEPERTURBATION

#include <array>
#include <vector>

namespace thames::perturbations::baseperturbation{
    
    ///////////
    // Reals //
    ///////////

    /**
     * @brief Class for the base perturbation.
     * 
     * @tparam T Numeric type.
     */
    template<class T>
    class BasePerturbation{

        private:

        public:

            /**
             * @brief Construct a new Base Perturbation object.
             * 
             */
            BasePerturbation();

            /**
             * @brief Destroy the Base Perturbation object.
             * 
             */
            ~BasePerturbation();

            ////////////
            // Arrays //
            ////////////

            /**
             * @brief Default total perturbing acceleration.
             * 
             * Returns zero total perturbing acceleration.
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::array<T, 3> Total perturbing acceleration.
             */
            virtual std::array<T, 3> acceleration_total(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const;

            /**
             * @brief Default non-potential perturbing acceleration.
             * 
             * Returns zero non-potential perturbing acceleration.
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::array<T, 3> Non-potential perturbing acceleration.
             */
            virtual std::array<T, 3> acceleration_nonpotential(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const;

            /**
             * @brief Default perturbing potential.
             * 
             * Returns zero potential.
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return T Perturbing potential.
             */
            virtual T potential(const T& t, const std::array<T, 3>& R) const;

            /**
             * @brief Default time derivative of the perturbing potential.
             * 
             * Returns zero time derivative of the perturbing potential.
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return T Time derivative of the perturbing potential.
             */
            virtual T potential_derivative(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const;

            /////////////
            // Vectors //
            /////////////

            /**
             * @brief Default total perturbing acceleration.
             * 
             * Returns zero total perturbing acceleration.
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
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     */
    template<class T, template<class> class P>
    class BasePerturbationPolynomial{

        private:

        public:

            /**
             * @brief Construct a new Base Perturbation Polynomial object.
             * 
             */
            BasePerturbationPolynomial();

            /**
             * @brief Destroy the Base Perturbation Polynomial object.
             * 
             */
            ~BasePerturbationPolynomial();

            /**
             * @brief Default total perturbing acceleration.
             * 
             * Returns zero total perturbing acceleration.
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