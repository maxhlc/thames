#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
using namespace smartuq::polynomial;
#endif

#include "cartesian.h"

namespace thames::conversions::cartesian {

    #ifdef THAMES_USE_SMARTUQ

    template<class T, template<class> class P>
    void cartesian_to_polynomial(const std::vector<T>& RV, const std::vector<T>& RVunc, int degree, std::vector<P<T>>& RVPolynomial) {
        // Clear polynomial vector
        RVPolynomial.clear();

        // Iterate through Cartesian state variables
        for(unsigned int ii=0; ii<6; ii++){
            // Generate polynomial
            RVPolynomial.push_back(P<T>(6, degree, ii, RV[ii]-RVunc[ii], RV[ii]+RVunc[ii]));

            // Initialise fast multiplication
            RVPolynomial[ii].initialize_M(6, degree);  
        }
          
    }
    template void cartesian_to_polynomial(const std::vector<double>& RV, const std::vector<double>& RVunc, int degree, std::vector<taylor_polynomial<double>>& RVPolynomial);
    template void cartesian_to_polynomial(const std::vector<double>& RV, const std::vector<double>& RVunc, int degree, std::vector<chebyshev_polynomial<double>>& RVPolynomial);

    #endif

}