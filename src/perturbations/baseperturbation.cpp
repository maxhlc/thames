#include <array>

#include "baseperturbation.h"

namespace thames::perturbations::baseperturbation{

    template<class T>
    BasePerturbation<T>::BasePerturbation(){

    };

    template<class T>
    BasePerturbation<T>::~BasePerturbation(){

    };
    
    template<class T>
    std::array<T, 3> BasePerturbation<T>::acceleration_total(T t, std::array<T, 3> R, std::array<T, 3> V) const{
        std::array<T, 3> F;
        F[0] = 0.0;
        F[1] = 0.0;
        F[2] = 0.0;
        return F;
    };

    template<class T>
    std::array<T, 3> BasePerturbation<T>::acceleration_nonpotential(T t, std::array<T, 3> R, std::array<T, 3> V) const{
        std::array<T, 3> F;
        F[0] = 0.0;
        F[1] = 0.0;
        F[2] = 0.0;
        return F;
    };

    template<class T>
    T BasePerturbation<T>::potential(T t, std::array<T, 3> R) const{
        return 0.0;
    }

    template<class T>
    T BasePerturbation<T>::potential_derivative(T t, std::array<T, 3> R, std::array<T, 3> V) const{
        return 0.0;
    }

    template class BasePerturbation<double>;

}