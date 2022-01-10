#ifndef THAMES_PERTURBATIONS_BASEPERTURBATION
#define THAMES_PERTURBATIONS_BASEPERTURBATION

#include "../types.h"

using namespace thames::types;

namespace thames::perturbations::baseperturbation{
    
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
            virtual Vector3 acceleration_total(double t, Vector3 R, Vector3 V);

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
            virtual Vector3 acceleration_nonpotential(double t, Vector3 R, Vector3 V);

            /**
             * @brief Default perturbing potential.
             * 
             * Returns zero potential.
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return double Perturbing potential.
             */
            virtual double potential(double t, Vector3 R);

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
            virtual double potential_derivative(double t, Vector3 R, Vector3 V);

    };

}

#endif