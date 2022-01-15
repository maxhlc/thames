#ifndef THAMES_UTIL_VECTOR
#define THAMES_UTIL_VECTOR

#include "../types.h"

using namespace thames::types;

namespace thames::util::vector{

    /**
     * @brief Function to calculate the dot product of two vectors with three elements.
     * 
     * @tparam real Type for real numbers (e.g. float, double, etc.)
     * @tparam vector Type for vector (e.g. std::vector<double, 3>, Eigen::Vector3d)
     * @param[in] a First vector.
     * @param[in] b Second vector.
     * @return real Dot product of the vectors.
     */
    template<class real, class vector>
    real dot3(const vector& a, const vector& b);

    /**
     * @brief Function to calculate the length of a vector with three elements.
     * 
     * @tparam real Type for real numbers (e.g. float, double, etc.)
     * @tparam vector Type for vector (e.g. std::vector<double, 3>, Eigen::Vector3d)
     * @param[in] a Vector.
     * @return real Length of the vector.
     */
    template<class real, class vector>
    real norm3(const vector& a);

    /**
     * @brief Function to return the slice of a vector.
     * 
     * @tparam vectorout Type for output vector.
     * @tparam vectorin Type for input vector.
     * @tparam integer Type for slicing bounds.
     * @param[in] v Input vector for slicing.
     * @param[in] a Lower index.
     * @param[in] b Upper index.
     * @return vectorout Output vector which has been sliced.
     */
    template<class vectorout, class vectorin, class integer>
    vectorout slice(const vectorin& v, const integer a, const integer b);

    template<class real, class vector>
    vector mult3(const real a, const vector& vec);

    template<class vector>
    vector cross3(const vector& a, const vector& b);

}

#endif