#ifndef THAMES_UTIL_POLYNOMIALS
#define THAMES_UTIL_POLYNOMIALS

#include <vector>

namespace thames::util::polynomials {

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Evaluate vector of polynomials at a point.
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] polynomials Vector of polynomials.
     * @param[in] x Evaluation point. 
     * @return std::vector<T> Point state vector.
     */
    template<class T, template<class> class P>
    std::vector<T> evaluate_polynomials(const std::vector<P<T>>& polynomials, const std::vector<T>& x);

    #endif

}

#endif