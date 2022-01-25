#include <array>
#include <vector>

#include "baseperturbation.h"

namespace thames::perturbations::baseperturbation{

    template<class T>
    BasePerturbation<T>::BasePerturbation(){

    };

    template<class T>
    BasePerturbation<T>::~BasePerturbation(){

    };

    ////////////
    // Arrays //
    ////////////
    
    template<class T>
    std::array<T, 3> BasePerturbation<T>::acceleration_total(T t, std::array<T, 3> R, std::array<T, 3> V) const{
        std::array<T, 3> F = {0.0, 0.0, 0.0};
        return F;
    };

    template<class T>
    std::array<T, 3> BasePerturbation<T>::acceleration_nonpotential(T t, std::array<T, 3> R, std::array<T, 3> V) const{
        std::array<T, 3> F = {0.0, 0.0, 0.0};
        return F;
    };

    template<class T>
    T BasePerturbation<T>::potential(T t, std::array<T, 3> R) const{
        T U = 0.0;
        return U;
    }

    template<class T>
    T BasePerturbation<T>::potential_derivative(T t, std::array<T, 3> R, std::array<T, 3> V) const{
        T Ut = 0.0;
        return Ut;
    }

    /////////////
    // Vectors //
    /////////////

    template<class T>
    std::vector<T> BasePerturbation<T>::acceleration_total(T t, std::vector<T> R, std::vector<T> V) const{
        std::vector<T> F = {0.0, 0.0, 0.0};
        return F;
    };

    template<class T>
    std::vector<T> BasePerturbation<T>::acceleration_nonpotential(T t, std::vector<T> R, std::vector<T> V) const{
        std::vector<T> F = {0.0, 0.0, 0.0};
        return F;
    };

    template<class T>
    T BasePerturbation<T>::potential(T t, std::vector<T> R) const{
        T U = 0.0;
        return U;
    }

    template<class T>
    T BasePerturbation<T>::potential_derivative(T t, std::vector<T> R, std::vector<T> V) const{
        T Ut = 0.0;
        return Ut;
    }

    template class BasePerturbation<double>;

}