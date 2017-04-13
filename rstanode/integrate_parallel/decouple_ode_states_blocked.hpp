#ifndef STAN_MATH_REV_ARR_FUN_DECOUPLE_ODE_STATES_BLOCKED_HPP
#define STAN_MATH_REV_ARR_FUN_DECOUPLE_ODE_STATES_BLOCKED_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <vector>

namespace stan {
  namespace math {

    /**
     * Takes sensitivity output from integrators and returns results
     * in precomputed_gradients format.
     *
     * Solution input vector size depends on requested sensitivities,
     * which can be enabled for initials and parameters. For each
     * sensitivity N states are computed. The input vector is expected
     * to be ordered by states, i.e. first the N states, then
     * optionally the N sensitivities for the initials (first the N
     * states for the first initial and so on), finally the
     * sensitivities for the M parameters follow optionally.
     *
     * @tparam T1_initial type of scalars for initial values.
     * @tparam T2_param type of scalars for parameters.
     * @param[in] y output from integrator in column-major format
     * as a coupled system output
     * @param[in] y0 initial state.
     * @param[in] theta parameters
     * @return a vector of states for each entry in y in Stan var
     * format
     */
    template <typename T_initial, typename T_param>
    inline
    std::vector<std::vector<typename stan::return_type<T_initial,
                                                       T_param>::type> >
    decouple_ode_states_blocked(const std::vector<std::vector<double> >& y,
                                const std::vector<std::vector<T_initial> >& y0,
                                const std::vector<std::vector<T_param> >& theta,
                                const std::vector<int>& Si_) {
      using std::vector;

      typedef stan::is_var<T_initial> initial_var;
      typedef stan::is_var<T_param> param_var;

      const size_t N = y0[0].size();
      const size_t M = theta[0].size();
      const size_t S = (initial_var::value ? N : 0)
        + (param_var::value ? M : 0);

      vector<var> temp_vars(N);
      vector<double> temp_gradients(S);
      vector<vector<var> > y_return(y.size());

      //size_t ic = 0;
      for (size_t s = 0; s < Si_.size()-1; ++s) {
        vector<typename stan::return_type<T_initial, T_param>::type> vars;
        vars.reserve(S);
        if (initial_var::value)
          vars.insert(vars.end(), y0[s].begin(), y0[s].end());
        if (param_var::value)
          vars.insert(vars.end(), theta[s].begin(), theta[s].end());

        for (size_t i = Si_[s]-1; i < Si_[s+1]-1; ++i) {
          for (size_t j = 0; j < N; ++j) {
            for (size_t k = 0; k < S; ++k) {
              //temp_gradients[k] = y[ic + i][N + N * k + j];
              temp_gradients[k] = y[i][N + N * k + j];
            }
            //temp_vars[j] = precomputed_gradients(y[ic + i][j],
            //                                     vars, temp_gradients);
            temp_vars[j] = precomputed_gradients(y[i][j],
                                                 vars, temp_gradients);
          }
          //y_return[ic + i] = temp_vars;
          y_return[i] = temp_vars;
        }

        //ic = Si_[s+1] - 1;
      }
      return y_return;
    }

    /**
     * The decouple ODE states operation for the case of no
     * sensitivities is equal to the indentity operation.
     *
     * @param[in] y output from integrator
     * @param[in] y0 initial state.
     * @param[in] theta parameters
     * @return y
     */
    template <>
    inline
    std::vector<std::vector<double> >
    decouple_ode_states_blocked(const std::vector<std::vector<double> >& y,
                                const std::vector<std::vector<double> >& y0,
                                const std::vector<std::vector<double> >& theta,
                                const std::vector<int>& Si_) {
      return y;
    }

  }
}
#endif
