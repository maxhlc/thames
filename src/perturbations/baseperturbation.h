#ifndef THAMES_PERTURBATIONS_BASEPERTURBATION
#define THAMES_PERTURBATIONS_BASEPERTURBATION

#include <array>
#include <vector>

namespace thames::perturbations::baseperturbation{
    
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
            virtual std::array<T, 3> acceleration_total(T t, std::array<T, 3> R, std::array<T, 3> V) const;

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
            virtual std::vector<T> acceleration_total(T t, std::vector<T> R, std::vector<T> V) const;

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
            virtual std::array<T, 3> acceleration_nonpotential(T t, std::array<T, 3> R, std::array<T, 3> V) const;

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
            virtual std::vector<T> acceleration_nonpotential(T t, std::vector<T> R, std::vector<T> V) const;

            /**
             * @brief Default perturbing potential.
             * 
             * Returns zero potential.
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return T Perturbing potential.
             */
            virtual T potential(T t, std::array<T, 3> R) const;

            /**
             * @brief Default perturbing potential.
             * 
             * Returns zero potential.
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return T Perturbing potential.
             */
            virtual T potential(T t, std::vector<T> R) const;

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
            virtual T potential_derivative(T t, std::array<T, 3> R, std::array<T, 3> V) const;

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
            virtual T potential_derivative(T t, std::vector<T> R, std::vector<T> V) const;

    };

}

#endif