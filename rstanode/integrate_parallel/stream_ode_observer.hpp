#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_STREAM_ODE_OBSERVER_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_STREAM_ODE_OBSERVER_HPP

#include <vector>

namespace stan {
  namespace math {

    /**
     * Observer for the coupled states.  Holds a reference to
     * an externally defined vector of vectors passed in at
     * construction time.
     */
    struct stream_ode_observer {
      std::vector<std::vector<double> >::iterator out_;
      int n_;

      /**
       * Construct a coupled ODE observer from the specified coupled
       * vector.
       *
       * @param y_coupled reference to a vector of vector of doubles.
       */
      explicit stream_ode_observer(std::vector<std::vector<double> >::iterator out)
        : out_(out), n_(0) {
      }

      /**
       * Callback function for Boost's ODE solver to record values.
       *
       * @param state solution at the specified time.
       * @param t time of solution.
       */
      void operator()(const std::vector<double>& state,
                      double t) {
        // we skip the first observation which is the initial value
        if(n_++ != 0) {
          *out_ = state;
          out_++;
        }
      }
    };

  }

}

#endif
