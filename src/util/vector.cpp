#include <cmath>

#include "vector.h"
#include "../types.h"

using namespace thames::types;

namespace thames::util::vector{

    template<class real, class vector>
    real dot3(const vector& a, const vector& b){
        // Return dot product
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    template double dot3<double, Vector3>(const Vector3&, const Vector3&);

    template<class real, class vector>
    real norm3(const vector& a){
        // Return square root of the dot product of the vector and itself
        return sqrt(dot3<real, vector>(a, a));
    }
    template double norm3<double, Vector3>(const Vector3&);

    template<class vectorout, class vectorin, class integer>
    vectorout slice(const vectorin& v, const integer a, const integer b){
        // Declare output vector
        vectorout vres;

        // Iterate through slice indices and populate output vector
        for(integer ii=a; ii<=b; ii++){
            vres[ii-a] = v[ii];
        }

        // Return output vector
        return vres;
    }
    template Vector3 slice<Vector3, Vector6, unsigned int>(const Vector6&, const unsigned int, const unsigned int);

    template<class real, class vector>
    vector mult3(const real a, const vector& vec){
        // Declare output vector
        vector vecout;

        // Interate through elements
        for(unsigned int ii=0; ii<3; ii++){
            vecout[ii] = a*vec[ii];
        }

        // Return output vector
        return vecout;
    }
    template Vector3 mult3<double, Vector3>(const double, const Vector3&);

    template<class vector>
    vector cross3(const vector& a, const vector& b){
        // Declare output vector
        vector vecout;

        // Calculate cross product
        vecout[0] = a[1]*b[2] - a[2]*b[1];
        vecout[1] = a[2]*b[0] - a[0]*b[2];
        vecout[2] = a[0]*b[1] - a[1]*b[0];

        // Return cross product
        return vecout;
    }
    template Vector3 cross3<Vector3>(const Vector3&, const Vector3&);

}