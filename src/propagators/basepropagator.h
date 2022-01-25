#ifndef THAMES_PROPAGATORS_BASEPROPAGATOR
#define THAMES_PROPAGATORS_BASEPROPAGATOR

#include "../perturbations/baseperturbation.h"

using namespace thames::perturbations::baseperturbation;

namespace thames::propagators::basepropagator {

    template<class T>
    class BasePropagator {

        private:

        public:

            ////////////
            // Arrays //
            ////////////

            virtual void derivative(const std::array<T, 6>& RV, std::array<T, 6>& RVdot, const T t) const;

            virtual std::array<T, 6> propagate(T tstart, T tend, T tstep, std::array<T, 6> RV, T atol, T rtol) const;

            /////////////
            // Vectors //
            /////////////

            virtual void derivative(const std::vector<T>& RV, std::vector<T>& RVdot, const T t) const;

            virtual std::vector<T> propagate(T tstart, T tend, T tstep, std::vector<T> RV, T atol, T rtol) const;           

    };

}

#endif