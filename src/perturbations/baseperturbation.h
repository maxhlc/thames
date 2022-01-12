#ifndef THAMES_PERTURBATIONS_BASEPERTURBATION
#define THAMES_PERTURBATIONS_BASEPERTURBATION

#include "../types.h"

using namespace thames::types;

namespace thames::perturbations::baseperturbation{
    
    /**
     * @brief Class for the base perturbation.
     * 
     * @tparam real Type for real numbers (e.g. float, double, etc.)
     * @tparam vector Type for vector (e.g. std::vector<double, 3>, Eigen::Vector3d)
     */
    template<class real, class vector>
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
             * @return vector Total perturbing acceleration.
             */
            virtual vector acceleration_total(real t, vector R, vector V);

            /**
             * @brief Default non-potential perturbing acceleration.
             * 
             * Returns zero non-potential perturbing acceleration.
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return vector Non-potential perturbing acceleration.
             */
            virtual vector acceleration_nonpotential(real t, vector R, vector V);

            /**
             * @brief Default perturbing potential.
             * 
             * Returns zero potential.
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return real Perturbing potential.
             */
            virtual real potential(real t, vector R);

            /**
             * @brief Default time derivative of the perturbing potential.
             * 
             * Returns zero time derivative of the perturbing potential.
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return real Time derivative of the perturbing potential.
             */
            virtual real potential_derivative(real t, vector R, vector V);

    };

}

#endif