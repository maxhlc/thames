#ifndef THAMES_CONVERSIONS_CARTESIAN
#define THAMES_CONVERSIONS_CARTESIAN

#include <vector>

namespace thames::conversions::cartesian {

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Generate vector of polynomials representing Cartesian state, and its uncertainty.
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] RV Cartesian state.
     * @param[in] RVunc Cartesian state uncertainty.
     * @param[in] degree Polynomial degree.
     * @param[out] RVPolynomial Polynomial Cartesian state.
     */
    template<class T, template<class> class P>
    void cartesian_to_polynomial(const std::vector<T>& RV, const std::vector<T>& RVunc, int degree, std::vector<P<T>>& RVPolynomial);

    #endif

}

#endif