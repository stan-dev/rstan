#ifndef STAN_MATH_REV_MAT_FUNCTOR_INTEGRATE_ODE_CVODES_STREAM_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_INTEGRATE_ODE_CVODES_STREAM_HPP

#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/err/check_consistent_size.hpp>
#include <stan/math/prim/arr/err/check_nonzero_size.hpp>
#include <stan/math/prim/arr/err/check_ordered.hpp>
#include <stan/math/rev/scal/meta/is_var.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/rev/mat/functor/cvodes_utils.hpp>
//#include <stan/math/rev/mat/functor/cvodes_ode_data.hpp>
#include <stan/math/rev/arr/fun/decouple_ode_states.hpp>
#include <cvodes/cvodes.h>
#include <cvodes/cvodes_band.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <algorithm>
#include <ostream>
#include <vector>

namespace stan {
  namespace math {

    /**
     * Free memory allocated for CVODES state, sensitivity, and
     * general memory.
     *
     * @param[in] cvodes_state State vector.
     * @param[in] cvodes_state_sens Sensivity vector.
     * @param[in] cvodes_mem Memory held for CVODES.
     * @param[in] S Number of sensitivities being calculated.
     */
    /*
    inline void free_cvodes_memory(N_Vector& cvodes_state,
                                   N_Vector* cvodes_state_sens,
                                   void* cvodes_mem, size_t S) {
      N_VDestroy_Serial(cvodes_state);
      if (cvodes_state_sens != NULL)
        N_VDestroyVectorArray_Serial(cvodes_state_sens, S);
      CVodeFree(&cvodes_mem);
    }
    */

    /**
     * Return the solutions for the specified system of ordinary
     * differential equations given the specified initial state,
     * initial times, times of desired solution, and parameters and
     * data, writing error and warning messages to the specified
     * stream.
     *
     * THIS VERSION IS openMP thread safe by not decoupling, but
     * returning the bare coupled system.
     *
     * This function is templated to allow the initial times to be
     * either data or autodiff variables and the parameters to be data
     * or autodiff variables.  The autodiff-based implementation for
     * reverse-mode are defined in namespace <code>stan::math</code>
     * and may be invoked via argument-dependent lookup by including
     * their headers.
     *
     * The solver used is based on the backward differentiation
     * formula which is an implicit numerical integration scheme
     * appropiate for stiff ODE systems.
     *
     * @tparam F type of ODE system function.
     * @tparam T_initial type of scalars for initial values.
     * @tparam T_param type of scalars for parameters.
     * @param[in] f functor for the base ordinary differential equation.
     * @param[in] y0 initial state.
     * @param[in] t0 initial time.
     * @param[in] ts times of the desired solutions, in strictly
     * increasing order, all greater than the initial time.
     * @param[in] theta parameter vector for the ODE.
     * @param[in] x continuous data vector for the ODE.
     * @param[in] x_int integer data vector for the ODE.
     * @param[in, out] msgs the print stream for warning messages.
     * @param[in] relative_tolerance relative tolerance passed to CVODE.
     * @param[in] absolute_tolerance absolute tolerance passed to CVODE.
     * @param[in] max_num_steps maximal number of admissable steps
     * between time-points
     * @return a vector of states, each state being a vector of the
     * same size as the state variable, corresponding to a time in ts.
     */
    template <typename F, typename T_initial, typename T_param>
    //std::vector<std::vector<typename stan::return_type<T_initial,
    //                                                   T_param>::type>
    //                                                   >
    void
    integrate_ode_cvodes_stream(const F& f,
                                const std::vector<T_initial>& y0,
                                const double t0,
                                const std::vector<double>& ts,
                                const std::vector<T_param>& theta,
                                const std::vector<double>& x,
                                const std::vector<int>& x_int,
                                const std::vector<double>& pbar, //scale
                                                                //information
                                                                //of
                                                                //sensitivity
                                                                //parameters
                                std::vector<std::vector<double> >::iterator out,
                                std::ostream* msgs = 0,
                                double relative_tolerance = 1e-10,
                                double absolute_tolerance = 1e-10,
                                long int max_num_steps = 1e8, // NOLINT(runtime/int)
                                int solver = 1) { // solver=1 stiff; =2 non-stiff
      typedef stan::is_var<T_initial> initial_var;
      typedef stan::is_var<T_param> param_var;

      check_finite("integrate_ode_cvodes", "initial state", y0);
      check_finite("integrate_ode_cvodes", "initial time", t0);
      check_finite("integrate_ode_cvodes", "times", ts);
      check_finite("integrate_ode_cvodes", "parameter vector", theta);
      check_finite("integrate_ode_cvodes", "continuous data", x);
      check_nonzero_size("integrate_ode_cvodes", "times", ts);
      check_nonzero_size("integrate_ode_cvodes", "initial state", y0);
      check_ordered("integrate_ode_cvodes", "times", ts);
      check_less("integrate_ode_cvodes", "initial time", t0, ts[0]);
      check_consistent_size("integrate_ode_cvodes", "pbar", pbar, y0.size() + theta.size());
      check_positive_finite("integrate_ode_cvodes", "pbar", pbar);
      if (relative_tolerance <= 0)
        invalid_argument("integrate_ode_cvodes",
                         "relative_tolerance,", relative_tolerance,
                         "", ", must be greater than 0");
      if (absolute_tolerance <= 0)
        invalid_argument("integrate_ode_cvodes",
                         "absolute_tolerance,", absolute_tolerance,
                         "", ", must be greater than 0");
      if (max_num_steps <= 0)
        invalid_argument("integrate_ode_cvodes",
                         "max_num_steps,", max_num_steps,
                         "", ", must be greater than 0");
      if (solver < 1 || solver > 2)
        invalid_argument("integrate_ode_cvodes",
                         "solver", solver,
                         "", ", must be 1/2 (BDF/Adams)");

                         const size_t N = y0.size();
      const size_t M = theta.size();
      // total number of sensitivities for initial values and params
      const size_t S = (initial_var::value ? N : 0)
        + (param_var::value ? M : 0);
      const size_t size = N * (S + 1);  // size of the coupled system
      std::vector<double> state(value_of(y0));
      N_Vector cvodes_state(N_VMake_Serial(N, &state[0]));
      N_Vector* cvodes_state_sens = NULL;

      // setup pbar vector acutally used
      std::vector<double> pbar_run;
      if (initial_var::value)
        pbar_run.insert(pbar_run.end(), pbar.begin(), pbar.begin() + N);
      if (param_var::value)
        pbar_run.insert(pbar_run.end(), pbar.begin() + N, pbar.end() );

      // PROBLEM!!!??? Maybe not as we are only passing references to
      // vars around and don't anything with them. BUT: Redesign such
      // that NO var ever is given to any ODE integrator.
      typedef cvodes_ode_data<F, T_initial, T_param> ode_data;
      ode_data cvodes_data(f, y0, theta, x, x_int, msgs);

      void* cvodes_mem = NULL;
      if(solver == 1)
        cvodes_mem = CVodeCreate(CV_BDF, CV_NEWTON);
      if(solver == 2)
        cvodes_mem = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
      if (cvodes_mem == NULL)
        throw std::runtime_error("CVodeCreate failed to allocate memory");

      //std::vector<std::vector<double> >
      //  y_coupled(ts.size(), std::vector<double>(size, 0));

      try {
        cvodes_check_flag(CVodeInit(cvodes_mem, &ode_data::ode_rhs,
                                    t0, cvodes_state),
                          "CVodeInit");

        // Assign pointer to this as user data
        cvodes_check_flag(CVodeSetUserData(cvodes_mem,
                            reinterpret_cast<void*>(&cvodes_data)),
                          "CVodeSetUserData");

        cvodes_set_options(cvodes_mem,
                           relative_tolerance, absolute_tolerance,
                           max_num_steps);

        // for the stiff solvers we need to reserve additional
        // memory and provide a Jacobian function call
        if(solver == 1) {
          cvodes_check_flag(CVDense(cvodes_mem, N), "CVDense");
          cvodes_check_flag(CVDlsSetDenseJacFn(cvodes_mem,
                                               &ode_data::dense_jacobian),
                            "CVDlsSetDenseJacFn");
        }

        // initialize forward sensitivity system of CVODES as needed
        if (S > 0) {
          cvodes_state_sens = N_VCloneVectorArray_Serial(S, cvodes_state);
          for (size_t s = 0; s < S; s++)
            N_VConst(RCONST(0.0), cvodes_state_sens[s]);

          // for varying initials, first N sensitivity systems
          // are for initials which have as initial the identity matrix
          if (initial_var::value) {
            for (size_t n = 0; n < N; n++)
              NV_Ith_S(cvodes_state_sens[n], n) = 1.0;
          }
          cvodes_check_flag(CVodeSensInit(cvodes_mem, static_cast<int>(S),
                                          CV_STAGGERED,
                                          &ode_data::ode_rhs_sens,
                                          cvodes_state_sens),
                            "CVodeSensInit");

          // provide scale information in parameters which deterime
          // absolute integration tolerances
          cvodes_check_flag(CVodeSetSensParams(cvodes_mem, NULL, &pbar_run[0], NULL),
                            "CVodeSetSensParams");

          cvodes_check_flag(CVodeSensEEtolerances(cvodes_mem),
                            "CVodeSensEEtolerances");
        }

        double t_init = t0;
        for (size_t n = 0; n < ts.size(); ++n) {
          double t_final = ts[n];
          std::vector<double> y_coupled(size);
          if (t_final != t_init)
            cvodes_check_flag(CVode(cvodes_mem, t_final, cvodes_state,
                                    &t_init, CV_NORMAL),
                              "CVode");
          std::copy(state.begin(), state.end(), y_coupled.begin());
          if (S > 0) {
            cvodes_check_flag(CVodeGetSens(cvodes_mem, &t_init,
                                           cvodes_state_sens),
                              "CVodeGetSens");
            for (size_t s = 0; s < S; s++)
              std::copy(NV_DATA_S(cvodes_state_sens[s]),
                        NV_DATA_S(cvodes_state_sens[s]) + N,
                        y_coupled.begin() + N + s * N);
          }
          *out = y_coupled;
          out++;
          t_init = t_final;
        }
      } catch (const std::exception& e) {
        free_cvodes_memory(cvodes_state, cvodes_state_sens, cvodes_mem, S);
        throw;
      }

      free_cvodes_memory(cvodes_state, cvodes_state_sens, cvodes_mem, S);

      // we return the bare sensitvity system; caller needs to
      // reformat to

      //std::vector<std::vector<typename stan::return_type<T_initial,
      //                                                   T_param>::type> > res;

      // the decoupling operation has to access the AD stack and hence
      // we need to instruct OpenMP to do things serially.
      //#pragma omp ordered
      //{
      //  res = decouple_ode_states(y_coupled, y0, theta);
      //}

      return;
    }

  }
}
#endif
