#include <Eigen/Core>

#include "vector.h"
#include "../types.h"

using namespace thames::types;

namespace thames::util::vector{

    template<class real, class vector>
    real dot3(vector a, vector b){
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    template double dot3<double, Vector3>(Vector3, Vector3);

    template<class real, class vector>
    real norm3(vector a){
        return sqrt(dot3<real, vector>(a, a));
    }
    template double norm3<double, Vector3>(Vector3);

    template<class vectorout, class vectorin, class integer>
    vectorout slice(vectorin v, integer a, integer b){
        vectorout vres;
        for(integer ii=a; ii<=b; ii++){
            vres[ii-a] = v[ii];
        }
        return vres;
    }
    template Vector3 slice<Vector3, Vector6, unsigned int>(Vector6, unsigned int, unsigned int);

}