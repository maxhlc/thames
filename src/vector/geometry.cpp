#include <array>
#include <cmath>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
using namespace smartuq::polynomial;
#endif

#include "geometry.h"

namespace thames::vector::geometry{

    ////////////
    // Arrays //
    ////////////

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

    /////////////
    // Vectors //
    /////////////

    template<class T>
    T dot3(const std::vector<T>& a, const std::vector<T>& b){
        // Return dot product
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    template double dot3<double>(const std::vector<double>&, const std::vector<double>&);

    template<class T>
    T norm3(const std::vector<T>& a){
        // Return square root of the dot product of the vector and itself
        return sqrt(dot3<T>(a, a));
    }
    template double norm3<double>(const std::vector<double>&);

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

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    template<class T, template<class> class P>
    P<T> dot3(const std::vector<P<T>>& a, const std::vector<P<T>>& b){
        // Return dot product
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    template taylor_polynomial<double> dot3(const std::vector<taylor_polynomial<double>>&, const std::vector<taylor_polynomial<double>>&);
    template chebyshev_polynomial<double> dot3(const std::vector<chebyshev_polynomial<double>>&, const std::vector<chebyshev_polynomial<double>>&);

    template<class T, template<class> class P>
    P<T> norm3(const std::vector<P<T>>& a){
        // Return square root of the dot product of the vector and itself
        return sqrt(dot3(a, a));
    }
    template taylor_polynomial<double> norm3(const std::vector<taylor_polynomial<double>>&);
    template chebyshev_polynomial<double> norm3(const std::vector<chebyshev_polynomial<double>>&);

    template<class T, template<class> class P>
    void cross3(const std::vector<P<T>>& a, const std::vector<P<T>>& b, std::vector<P<T>>& vecout){
        // Calculate cross product
        vecout[0] = a[1]*b[2] - a[2]*b[1];
        vecout[1] = a[2]*b[0] - a[0]*b[2];
        vecout[2] = a[0]*b[1] - a[1]*b[0];
    }
    template void cross3(const std::vector<taylor_polynomial<double>>&, const std::vector<taylor_polynomial<double>>&, std::vector<taylor_polynomial<double>>&);
    template void cross3(const std::vector<chebyshev_polynomial<double>>&, const std::vector<chebyshev_polynomial<double>>&, std::vector<chebyshev_polynomial<double>>&);

    template<class T, template<class> class P>
    std::vector<P<T>> cross3(const std::vector<P<T>>& a, const std::vector<P<T>>& b){
        // Declare output vector
        std::vector<P<T>> vecout(a);

        // Calculate cross product
        cross3(a, b, vecout);

        // Return cross product
        return vecout;
    }
    template std::vector<taylor_polynomial<double>> cross3(const std::vector<taylor_polynomial<double>>&, const std::vector<taylor_polynomial<double>>&);
    template std::vector<chebyshev_polynomial<double>> cross3(const std::vector<chebyshev_polynomial<double>>&, const std::vector<chebyshev_polynomial<double>>&);

    #endif

}