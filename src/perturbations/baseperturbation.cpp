#include "baseperturbation.h"
#include "../types.h"

using namespace thames::types;

namespace thames::perturbations::baseperturbation{

    template<class real, class vector>
    BasePerturbation<real, vector>::BasePerturbation(){

    };

    template<class real, class vector>
    BasePerturbation<real, vector>::~BasePerturbation(){

    };
    
    template<class real, class vector>
    vector BasePerturbation<real, vector>::acceleration_total(real t, vector R, vector V){
        vector F;
        F[0] = 0.0;
        F[1] = 0.0;
        F[2] = 0.0;
        return F;
    };

    template<class real, class vector>
    vector BasePerturbation<real, vector>::acceleration_nonpotential(real t, vector R, vector V){
        vector F;
        F[0] = 0.0;
        F[1] = 0.0;
        F[2] = 0.0;
        return F;
    };

    template<class real, class vector>
    real BasePerturbation<real, vector>::potential(real t, vector R){
        return 0.0;
    }

    template<class real, class vector>
    real BasePerturbation<real, vector>::potential_derivative(real t, vector R, vector V){
        return 0.0;
    }

    template class BasePerturbation<double, Vector3>;

}