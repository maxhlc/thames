#include <array>
#include <cmath>

#include "geometry.h"

namespace thames::vector::geometry{

    template<class T>
    T dot3(const std::array<T, 3>& a, const std::array<T, 3>& b){
        // Return dot product
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    template double dot3<double>(const std::array<double, 3>&, const std::array<double, 3>&);

    template<class T>
    T dot3(const std::vector<T>& a, const std::vector<T>& b){
        // Return dot product
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    template double dot3<double>(const std::vector<double>&, const std::vector<double>&);

    template<class T>
    T norm3(const std::array<T, 3>& a){
        // Return square root of the dot product of the vector and itself
        return sqrt(dot3<T>(a, a));
    }
    template double norm3<double>(const std::array<double, 3>&);

    template<class T>
    T norm3(const std::vector<T>& a){
        // Return square root of the dot product of the vector and itself
        return sqrt(dot3<T>(a, a));
    }
    template double norm3<double>(const std::vector<double>&);

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
    std::vector<T> slice(const std::vector<T>& v, const unsigned int a, const unsigned int b){
        // Declare output vector
        std::vector<T> vres(b-a+1);

        // Iterate through slice indices and populate output vector
        for(unsigned int ii=a; ii<=b; ii++){
            vres[ii-a] = v[ii];
        }

        // Return output vector
        return vres;
    }
    template std::vector<double> slice<double>(const std::vector<double>&, const unsigned int, const unsigned int);

    template<class T>
    void cross3(const std::array<T, 3>& a, const std::array<T, 3>& b, std::array<T, 3>& vecout){
        // Calculate cross product
        vecout[0] = a[1]*b[2] - a[2]*b[1];
        vecout[1] = a[2]*b[0] - a[0]*b[2];
        vecout[2] = a[0]*b[1] - a[1]*b[0];
    }
    template void cross3<double>(const std::array<double, 3>&, const std::array<double, 3>&, std::array<double, 3>&);

    template<class T>
    std::array<T, 3> cross3(const std::array<T, 3>& a, const std::array<T, 3>& b){
        // Declare output vector
        std::array<T, 3> vecout;

        // Calculate cross product
        cross3<T>(a, b, vecout);

        // Return cross product
        return vecout;
    }
    template std::array<double, 3> cross3<double>(const std::array<double, 3>&, const std::array<double, 3>&);

    template<class T>
    void cross3(const std::vector<T>& a, const std::vector<T>& b, std::vector<T>& vecout){
        // Calculate cross product
        vecout[0] = a[1]*b[2] - a[2]*b[1];
        vecout[1] = a[2]*b[0] - a[0]*b[2];
        vecout[2] = a[0]*b[1] - a[1]*b[0];
    }
    template void cross3<double>(const std::vector<double>&, const std::vector<double>&, std::vector<double>&);

    template<class T>
    std::vector<T> cross3(const std::vector<T>& a, const std::vector<T>& b){
        // Declare output vector
        std::vector<T> vecout(3);

        // Calculate cross product
        cross3<T>(a, b, vecout);

        // Return cross product
        return vecout;
    }
    template std::vector<double> cross3<double>(const std::vector<double>&, const std::vector<double>&);

}