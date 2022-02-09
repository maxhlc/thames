#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
using namespace smartuq::polynomial;
#endif

#include "polynomials.h"

namespace thames::util::polynomials {

    #ifdef THAMES_USE_SMARTUQ

    template<class T, template<class> class P>
    std::vector<T> evaluate_polynomials(const std::vector<P<T>>& polynomials, const std::vector<T>& x) {
        // Declare point state vector
        std::vector<T> numeric;

        // Evaluate each polynomial
        for(std::size_t ii=0; ii<polynomials.size(); ii++)
            numeric.push_back(polynomials[ii].evaluate(x));

        // Return point state vector
        return numeric;
    }
    template std::vector<double> evaluate_polynomials(const std::vector<taylor_polynomial<double>>&, const std::vector<double>&);
    template std::vector<double> evaluate_polynomials(const std::vector<chebyshev_polynomial<double>>&, const std::vector<double>&);

    #endif

}