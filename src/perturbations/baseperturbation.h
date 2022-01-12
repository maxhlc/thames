#ifndef THAMES_PERTURBATIONS_BASEPERTURBATION
#define THAMES_PERTURBATIONS_BASEPERTURBATION

#include "../types.h"

using namespace thames::types;

namespace thames::perturbations::baseperturbation{
    
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
             * @return Vector3 Total perturbing acceleration.
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
             * @return Vector3 Non-potential perturbing acceleration.
             */
            virtual vector acceleration_nonpotential(real t, vector R, vector V);

            /**
             * @brief Default perturbing potential.
             * 
             * Returns zero potential.
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return double Perturbing potential.
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
             * @return double Time derivative of the perturbing potential.
             */
            virtual real potential_derivative(real t, vector R, vector V);

    };

}

#endif