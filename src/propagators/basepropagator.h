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
             * @param[in] x State.
             * @param[out] dxdt State derivative.
             * @param[in] t Current time.
             */
            virtual void derivative(const std::array<T, 6>& x, std::array<T, 6>& dxdt, const T t) const;

            /**
             * @brief Propagation method.
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
             * @param[in] x State.
             * @param[out] dxdt State derivative.
             * @param[in] t Time.
             */
            virtual void derivative(const std::vector<T>& x, std::vector<T>& dxdt, const T t) const;

            /**
             * @brief Propagation method.
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
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     */
    template<class T, template<class> class P>
    class BasePropagatorPolynomial {

        private:

        public:

            // /**
            //  * @brief State derivative method.
            //  * 
            //  * @param[in] x State.
            //  * @param[out] dxdt State derivative.
            //  * @param[in] t Time.
            //  */
            // virtual void derivative(const std::vector<P<T>>& x, std::vector<P<T>>& dxdt, const T t) const;

            /**
             * @brief Propagation method.
             * 
             * @param[in] tstart Propagation start time in physical time.
             * @param[in] tend Propagation end time in physical time.
             * @param[in] tstep Initial timestep for propagation.
             * @param[in] RV Initial Cartesian state.
             * @return std::vector<P<T>> Final state.
             */
            virtual std::vector<P<T>> propagate(T tstart, T tend, T tstep, std::vector<P<T>> RV) const;           

    };

    #endif

}

#endif