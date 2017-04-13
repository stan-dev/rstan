// defines the integrate_parallel function which uses OpenMP to
// evaluate the ODE in parallel

namespace MODEL_NAMESPACE {

  template <typename T0, typename T1, typename T2, typename T4, typename T5, typename T9, typename T10, typename T13>
  std::vector<std::vector<typename boost::math::tools::promote_args<T0, T1, T2, T4, typename boost::math::tools::promote_args<T5, T9, T10, T13>::type>::type> >
  integrate_parallel(const std::vector<std::vector<T0> >& y0,
                     const std::vector<T1>& t0,
                     const std::vector<T2>& ts,
                     const std::vector<int>& Si_ts,
                     const std::vector<std::vector<T4> >& theta,
                     const std::vector<T5>& x_r, // NOTE:
                     // x_r must
                     // be given
                     // as data!
                     const std::vector<int>& Si_x_r,
                     const std::vector<int>& x_i,
                     const std::vector<int>& Si_x_i,
                     const T9& rel_tol, const T10& abs_tol, const int& max_steps,
                     const int& solver,
                     const std::vector<T13>& pbar,
                     std::ostream* pstream__) {
    int O = y0.size();
    int exceptions = 0;
    std::vector<std::vector<double> > y_coupled(ts.size());

#pragma omp parallel shared(exceptions)
    {
      //std::cout << "Hello from thread " << omp_get_thread_num() <<
      // ", nthreads " << omp_get_num_threads() << std::endl;

#pragma omp for schedule(runtime)
      for(int o = 0; o < O; ++o) {
        //std::vector<std::vector<double> > y_coupled_run;
        //y_coupled_run.reserve(S_ts[o]);
        try {
          if(exceptions == 0) {
            if(solver == 0) {
              stan::math::integrate_ode_rk45_stream(ode_functor(),
                                                    y0[o], t0[o],
                                                    stan::model::rvalue(ts, stan::model::cons_list(stan::model::index_min_max(Si_ts[o], Si_ts[o+1]-1), stan::model::nil_index_list()), "ts"),
                                                    theta[o],
                                                    stan::model::rvalue(x_r, stan::model::cons_list(stan::model::index_min_max(Si_x_r[o], Si_x_r[o+1]-1), stan::model::nil_index_list()), "x_r"),
                                                    stan::model::rvalue(x_i, stan::model::cons_list(stan::model::index_min_max(Si_x_i[o], Si_x_i[o+1]-1), stan::model::nil_index_list()), "x_i"),
                                                    y_coupled.begin() + Si_ts[o]-1,
                                                    pstream__,
                                                    rel_tol, abs_tol, max_steps);
            } else {
              stan::math::integrate_ode_cvodes_stream(ode_functor(),
                                                      y0[o], t0[o],
                                                      stan::model::rvalue(ts, stan::model::cons_list(stan::model::index_min_max(Si_ts[o], Si_ts[o+1]-1), stan::model::nil_index_list()), "ts"),
                                                      theta[o],
                                                      stan::model::rvalue(x_r, stan::model::cons_list(stan::model::index_min_max(Si_x_r[o], Si_x_r[o+1]-1), stan::model::nil_index_list()), "x_r"),
                                                      stan::model::rvalue(x_i, stan::model::cons_list(stan::model::index_min_max(Si_x_i[o], Si_x_i[o+1]-1), stan::model::nil_index_list()), "x_i"),
                                                      pbar,
                                                      y_coupled.begin() + Si_ts[o]-1,
                                                      pstream__,
                                                      rel_tol, abs_tol, max_steps,
                                                      solver);
            }
          }
        } catch(const std::exception& e) {
#pragma omp atomic
          ++exceptions;
        }
      }
    }
    if (exceptions != 0)
      throw std::domain_error("ODE error");

    return(stan::math::decouple_ode_states_blocked(y_coupled, y0, theta, Si_ts));

  }

}
