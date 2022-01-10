#include "baseperturbation.h"
#include "../types.h"

using namespace thames::types;

namespace thames::perturbations::baseperturbation{

    BasePerturbation::BasePerturbation(){

    };

    BasePerturbation::~BasePerturbation(){

    };
    
    Vector3 BasePerturbation::acceleration_total(double t, Vector3 R, Vector3 V){
        Vector3 F;
        F << 0.0, 0.0, 0.0;
        return F;
    };

    Vector3 BasePerturbation::acceleration_nonpotential(double t, Vector3 R, Vector3 V){
        Vector3 F;
        F << 0.0, 0.0, 0.0;
        return F;
    };

    double BasePerturbation::potential(double t, Vector3 R){
        return 0.0;
    }

    double BasePerturbation::potential_derivative(double t, Vector3 R, Vector3 V){
        return 0.0;
    }

}