/*
MIT License

Copyright (c) 2021-2022 Max Hallgarten La Casta

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <array>
#include <memory>
#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../../include/constants/statetypes.h"
#include "../../include/conversions/dimensional.h"
#include "../../include/conversions/geqoe.h"
#include "../../include/conversions/keplerian.h"
#include "../../include/conversions/universal.h"
#include "../../include/perturbations/baseperturbation.h"

namespace thames::conversions::universal {

    using thames::constants::statetypes::StateTypes;
    using thames::constants::statetypes::CARTESIAN;
    using thames::constants::statetypes::GEQOE;
    using thames::constants::statetypes::KEPLERIAN;
    using thames::conversions::dimensional::DimensionalFactors;
    using thames::perturbations::baseperturbation::BasePerturbation;

    ///////////
    // Reals //
    ///////////

    template<class T>
    std::vector<T> convert_state(const T& t, const std::vector<T>& state, const T& mu, const StateTypes& statetype1, const StateTypes& statetype2, const std::shared_ptr<const BasePerturbation<T>> perturbation) {
        // Return input directly if the two types match
        if (statetype1 == statetype2)
            return state;

        // Cartesian -> GEqOE
        if (statetype1 == CARTESIAN && statetype2 == GEQOE)
            return thames::conversions::geqoe::cartesian_to_geqoe(t, state, mu, perturbation);

        // GEqOE -> Cartesian
        if (statetype1 == GEQOE && statetype2 == CARTESIAN)
            return thames::conversions::geqoe::geqoe_to_cartesian(t, state, mu, perturbation);

        // Cartesian -> Keplerian
        if (statetype1 == CARTESIAN && statetype2 == KEPLERIAN)
            return thames::conversions::keplerian::cartesian_to_keplerian(state, mu);

        // Keplerian -> Cartesian
        if (statetype1 == KEPLERIAN && statetype2 == CARTESIAN)
            return thames::conversions::keplerian::keplerian_to_cartesian(state, mu);

        // Throw error if combinations not accepted
        throw std::runtime_error("Unsupported state conversion");
    }
    template std::vector<double> convert_state(const double&, const std::vector<double>&, const double&, const StateTypes&, const StateTypes&, const std::shared_ptr<const BasePerturbation<double>> perturbation);

    template<class T>
    std::vector<T> nondimensionalise_state(const std::vector<T>& state, const StateTypes& statetype, const DimensionalFactors<T>& factors) {
        // Switch through state types
        switch (statetype) {
            case CARTESIAN:
                return thames::conversions::dimensional::cartesian_nondimensionalise(state, factors);
                break;
            
            case GEQOE:
                return thames::conversions::dimensional::geqoe_nondimensionalise(state, factors);
                break;
        
            default:
                throw std::runtime_error("Unsupported non-dimensionalisation");
                break;
        }
    }
    template std::vector<double> nondimensionalise_state(const std::vector<double>& state, const StateTypes& statetype, const DimensionalFactors<double>& factors);

    template<class T>
    std::vector<T> dimensionalise_state(const std::vector<T>& statend, const StateTypes& statetype, const DimensionalFactors<T>& factors) {
        // Switch through state types
        switch (statetype) {
            case CARTESIAN:
                return thames::conversions::dimensional::cartesian_dimensionalise(statend, factors);
                break;
            
            case GEQOE:
                return thames::conversions::dimensional::geqoe_dimensionalise(statend, factors);
                break;
        
            default:
                throw std::runtime_error("Unsupported dimensionalisation");
                break;
        }
    }
    template std::vector<double> dimensionalise_state(const std::vector<double>& statend, const StateTypes& statetype, const DimensionalFactors<double>& factors);

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::polynomial;

    using thames::perturbations::baseperturbation::BasePerturbationPolynomial;

    template<class T, template <class> class P>
    std::vector<P<T>> convert_state(const T& t, const std::vector<P<T>>& state, const T& mu, const StateTypes& statetype1, const StateTypes& statetype2, const std::shared_ptr<const BasePerturbationPolynomial<T, P>> perturbation) {
        // Return input directly if the two types match
        if (statetype1 == statetype2)
            return state;

        // Cartesian -> GEqOE
        if (statetype1 == CARTESIAN && statetype2 == GEQOE)
            return thames::conversions::geqoe::cartesian_to_geqoe(t, state, mu, perturbation);

        // GEqOE -> Cartesian
        if (statetype1 == GEQOE && statetype2 == CARTESIAN)
            return thames::conversions::geqoe::geqoe_to_cartesian(t, state, mu, perturbation);

        // Throw error if combinations not accepted
        throw std::runtime_error("Unsupported state conversion");
    }
    template std::vector<taylor_polynomial<double>> convert_state(const double&, const std::vector<taylor_polynomial<double>>&, const double&, const StateTypes&, const StateTypes&, const std::shared_ptr<const BasePerturbationPolynomial<double, taylor_polynomial>>);
    template std::vector<chebyshev_polynomial<double>> convert_state(const double&, const std::vector<chebyshev_polynomial<double>>&, const double&, const StateTypes&, const StateTypes&, const std::shared_ptr<const BasePerturbationPolynomial<double, chebyshev_polynomial>>);

    template<class T, template <class> class P>
    std::vector<P<T>> nondimensionalise_state(const std::vector<P<T>>& state, const StateTypes& statetype, const DimensionalFactors<T>& factors) {
        // Switch through state types
        switch (statetype) {
            case CARTESIAN:
                return thames::conversions::dimensional::cartesian_nondimensionalise(state, factors);
                break;
            
            case GEQOE:
                return thames::conversions::dimensional::geqoe_nondimensionalise(state, factors);
                break;
        
            default:
                throw std::runtime_error("Unsupported non-dimensionalisation");
                break;
        }
    }
    template std::vector<taylor_polynomial<double>> nondimensionalise_state(const std::vector<taylor_polynomial<double>>& state, const StateTypes& statetype, const DimensionalFactors<double>& factors);
    template std::vector<chebyshev_polynomial<double>> nondimensionalise_state(const std::vector<chebyshev_polynomial<double>>& state, const StateTypes& statetype, const DimensionalFactors<double>& factors);

    template<class T, template <class> class P>
    std::vector<P<T>> dimensionalise_state(const std::vector<P<T>>& statend, const StateTypes& statetype, const DimensionalFactors<T>& factors) {
        // Switch through state types
        switch (statetype) {
            case CARTESIAN:
                return thames::conversions::dimensional::cartesian_dimensionalise(statend, factors);
                break;
            
            case GEQOE:
                return thames::conversions::dimensional::geqoe_dimensionalise(statend, factors);
                break;
        
            default:
                throw std::runtime_error("Unsupported dimensionalisation");
                break;
        }
    }
    template std::vector<taylor_polynomial<double>> dimensionalise_state(const std::vector<taylor_polynomial<double>>& statend, const StateTypes& statetype, const DimensionalFactors<double>& factors);
    template std::vector<chebyshev_polynomial<double>> dimensionalise_state(const std::vector<chebyshev_polynomial<double>>& statend, const StateTypes& statetype, const DimensionalFactors<double>& factors);

    #endif

}