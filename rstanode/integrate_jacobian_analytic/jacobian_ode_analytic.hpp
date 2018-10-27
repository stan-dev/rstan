
namespace stan {
namespace math {

template <>
void jacobian_ode_states_parameters(const ode_functor& f, const std::vector<double>& y,
                                    double t, const std::vector<double>& theta_,
                                    const std::vector<double>& x,
                                    const std::vector<int>& x_int,
                                    std::ostream* msgs, vector_d& dydt,
                                    matrix_d& Jy, matrix_d& Jtheta) {
  using std::pow;
  using std::exp;
        
#include ODE_HEADER(defs)
#include ODE_HEADER(rhs)
#include ODE_HEADER(jac_states)
#include ODE_HEADER(jac_params)
  
}

template <>
void jacobian_ode_states(const ode_functor& f, const std::vector<double>& y,
                         double t, const std::vector<double>& theta_,
                         const std::vector<double>& x,
                         const std::vector<int>& x_int,
                         std::ostream* msgs, vector_d& dydt,
                         matrix_d& Jy) {
  using std::pow;
  using std::exp;
        
#include ODE_HEADER(defs)
#include ODE_HEADER(rhs)
#include ODE_HEADER(jac_states)
  
}

}
}
