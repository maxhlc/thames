#include <array>
#include <cmath>

#include "vector.h"

namespace thames::util::vector{

    template<class T>
    T dot3(const std::array<T, 3>& a, const std::array<T, 3>& b){
        // Return dot product
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    template double dot3<double>(const std::array<double, 3>&, const std::array<double, 3>&);

    template<class T>
    T norm3(const std::array<T, 3>& a){
        // Return square root of the dot product of the vector and itself
        return sqrt(dot3<T>(a, a));
    }
    template double norm3<double>(const std::array<double, 3>&);

    template<class T, const std::size_t Si, const std::size_t So>
    std::array<T, So> slice(const std::array<T, Si>& v, const unsigned int a, const unsigned int b){
        // Declare output vector
        std::array<T, So> vres;

        // Iterate through slice indices and populate output vector
        for(unsigned int ii=a; ii<=b; ii++){
            vres[ii-a] = v[ii];
        }

        // Return output vector
        return vres;
    }
    template std::array<double, 3> slice<double, 6, 3>(const std::array<double, 6>&, const unsigned int, const unsigned int);

    template<class T>
    std::array<T, 3> mult3(const T a, const std::array<T, 3>& vec){
        // Declare output vector
        std::array<T, 3> vecout;

        // Interate through elements
        for(unsigned int ii=0; ii<3; ii++){
            vecout[ii] = a*vec[ii];
        }

        // Return output vector
        return vecout;
    }
    template std::array<double, 3> mult3<double>(const double, const std::array<double, 3>&);

    template<class T>
    std::array<T, 3> cross3(const std::array<T, 3>& a, const std::array<T, 3>& b){
        // Declare output vector
        std::array<T, 3> vecout;

        // Calculate cross product
        vecout[0] = a[1]*b[2] - a[2]*b[1];
        vecout[1] = a[2]*b[0] - a[0]*b[2];
        vecout[2] = a[0]*b[1] - a[1]*b[0];

        // Return cross product
        return vecout;
    }
    template std::array<double, 3> cross3<double>(const std::array<double, 3>&, const std::array<double, 3>&);

}