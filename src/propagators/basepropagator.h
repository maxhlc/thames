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

#ifndef THAMES_PROPAGATORS_BASEPROPAGATOR
#define THAMES_PROPAGATORS_BASEPROPAGATOR

#include <array>
#include <vector>

namespace thames::propagators::basepropagator {

    ///////////
    // Reals //
    ///////////

    /**
     * @brief Base propagator abstract object.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-01-26
     * 
     * @tparam T Numeric type.
     */
    template<class T>
    class BasePropagator {

        private:

        public:

            ////////////
            // Arrays //
            ////////////

            /**
             * @brief State derivative method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-26
             * 
             * @param[in] x State.
             * @param[out] dxdt State derivative.
             * @param[in] t Current time.
             */
            virtual void derivative(const std::array<T, 6>& x, std::array<T, 6>& dxdt, const T t) const;

            /**
             * @brief Propagation method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-26
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] RV Initial Cartesian state.
             * @param[in] atol Solver absolute tolerance.
             * @param[in] rtol Solver relative tolerance.
             * @return std::array<T, 6> Final state.
             */
            virtual std::array<T, 6> propagate(T tstart, T tend, T tstep, std::array<T, 6> RV, T atol, T rtol) const;

            /////////////
            // Vectors //
            /////////////

            /**
             * @brief State derivative method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-26
             * 
             * @param[in] x State.
             * @param[out] dxdt State derivative.
             * @param[in] t Time.
             */
            virtual void derivative(const std::vector<T>& x, std::vector<T>& dxdt, const T t) const;

            /**
             * @brief Propagation method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-01-26
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] RV Initial Cartesian state.
             * @param[in] atol Solver absolute tolerance.
             * @param[in] rtol Solver relative tolerance.
             * @return std::vector<T> Final state.
             */
            virtual std::vector<T> propagate(T tstart, T tend, T tstep, std::vector<T> RV, T atol, T rtol) const;           

    };

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Base propagator abstract object for polynomial propagations.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-02-22
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     */
    template<class T, template<class> class P>
    class BasePropagatorPolynomial {

        private:

        public:

            /**
             * @brief Propagation method.
             * 
             * @author Max Hallgarten La Casta
             * @date 2022-02-22
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] RV Initial Cartesian state.
             * @param[in] atol Solver absolute tolerance.
             * @param[in] rtol Solver relative tolerance.
             * @return std::vector<P<T>> Final state.
             */
            virtual std::vector<P<T>> propagate(T tstart, T tend, T tstep, std::vector<P<T>> RV, T atol, T rtol) const;         

    };

    #endif

}

#endif