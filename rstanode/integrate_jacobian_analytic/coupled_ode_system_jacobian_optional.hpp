#ifndef STAN_MATH_REV_ARR_FUNCTOR_COUPLED_ODE_SYSTEM_JACOBIAN_OPTIONAL_HPP
#define STAN_MATH_REV_ARR_FUNCTOR_COUPLED_ODE_SYSTEM_JACOBIAN_OPTIONAL_HPP

#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/arr/meta/get.hpp>
#include <stan/math/prim/arr/meta/length.hpp>
#include <stan/math/prim/arr/functor/coupled_ode_system.hpp>
#include <stan/math/prim/arr/fun/value_of.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/core.hpp>

#include <ostream>
#include <stdexcept>
#include <vector>

namespace stan {
namespace math {

template <typename F>
void jacobian_ode_states_parameters(const F& f, const std::vector<double>& z,
                                    double t, const std::vector<double>& theta,
                                    const std::vector<double>& x,
                                    const std::vector<int>& x_int,
                                    std::ostream* msgs,
                                    Eigen::Map<vector_d>& dy_dt,
                                    Eigen::Map<matrix_d>& Jy,
                                    Eigen::Map<matrix_d>& Jtheta) {
  using std::vector;
  const size_t N = dy_dt.size();
  const size_t M = theta.size();

  try {
    start_nested();

    vector<double> grad(N + M);
    Eigen::Map<vector_d> grad_y(grad.data(), N);
    Eigen::Map<vector_d> grad_theta(grad.data() + N, M);

    vector<var> z_vars;
    z_vars.reserve(N + M);

    vector<var> y_vars(z.begin(), z.begin() + N);
    z_vars.insert(z_vars.end(), y_vars.begin(), y_vars.end());

    vector<var> theta_vars(theta.begin(), theta.end());
    z_vars.insert(z_vars.end(), theta_vars.begin(), theta_vars.end());

    vector<var> dy_dt_vars = f(t, y_vars, theta_vars, x, x_int, msgs);

    check_size_match("coupled_ode_system", "dz_dt", dy_dt_vars.size(), "states",
                     N);

    for (size_t i = 0; i < N; i++) {
      dy_dt(i) = dy_dt_vars[i].val();
      set_zero_all_adjoints_nested();
      dy_dt_vars[i].grad(z_vars, grad);

      Jy.row(i) = grad_y;
      Jtheta.row(i) = grad_theta;
    }
  } catch (const std::exception& e) {
    recover_memory_nested();
    throw;
  }
  recover_memory_nested();
}

template <typename F>
void jacobian_ode_states(const F& f, const std::vector<double>& z, double t,
                         const std::vector<double>& theta,
                         const std::vector<double>& x,
                         const std::vector<int>& x_int,
                         std::ostream* msgs,
                         Eigen::Map<vector_d>& dy_dt,
                         Eigen::Map<matrix_d>& Jy
                         ) {
  using std::vector;
  const size_t N = dy_dt.size();
  //const size_t M = theta.size();

  try {
    start_nested();

    vector<double> grad(N);
    Eigen::Map<vector_d> grad_y(grad.data(), N);

    vector<var> y_vars(z.begin(), z.begin() + N);

    vector<var> dy_dt_vars = f(t, y_vars, theta, x, x_int, msgs);

    check_size_match("coupled_ode_system", "dz_dt", dy_dt_vars.size(), "states",
                     N);

    for (size_t i = 0; i < N; i++) {
      dy_dt(i) = dy_dt_vars[i].val();
      set_zero_all_adjoints_nested();
      dy_dt_vars[i].grad(y_vars, grad);

      Jy.row(i) = grad_y;
    }
  } catch (const std::exception& e) {
    recover_memory_nested();
    throw;
  }
  recover_memory_nested();
}

// This code is in this directory because it includes var
// It is in namespace stan::math so that the partial template
// specializations are treated as such.

/**
 * The coupled ODE system for known initial values and unknown
 * parameters.
 *
 * <p>If the base ODE state is size N and there are M parameters,
 * the coupled system has N + N * M states.
 * <p>The first N states correspond to the base system's N states:
 * \f$ \frac{d x_n}{dt} \f$
 *
 * <p>The next M states correspond to the sensitivities of the
 * parameters with respect to the first base system equation:
 * \f[
 *   \frac{d x_{N+m}}{dt}
 *   = \frac{d}{dt} \frac{\partial x_1}{\partial \theta_m}
 * \f]
 *
 * <p>The final M states correspond to the sensitivities with respect
 * to the second base system equation, etc.
 *
 * @tparam F type of functor for the base ode system.
 */
template <typename F>
struct coupled_ode_system<F, double, var> {
  const F& f_;
  const std::vector<double>& y0_dbl_;
  const std::vector<var>& theta_;
  const std::vector<double> theta_dbl_;
  // mutable std::vector<fvar<double>> theta_fvars_;
  const std::vector<double>& x_;
  const std::vector<int>& x_int_;
  const size_t N_;
  const size_t M_;
  const size_t size_;
  std::ostream* msgs_;

  /**
   * Construct a coupled ODE system with the specified base
   * ODE system, base initial state, parameters, data, and a
   * message stream.
   *
   * @param[in] f the base ODE system functor.
   * @param[in] y0 the initial state of the base ode.
   * @param[in] theta parameters of the base ode.
   * @param[in] x real data.
   * @param[in] x_int integer data.
   * @param[in, out] msgs stream to which messages are printed.
   */
  coupled_ode_system(const F& f, const std::vector<double>& y0,
                     const std::vector<var>& theta,
                     const std::vector<double>& x,
                     const std::vector<int>& x_int, std::ostream* msgs)
      : f_(f),
        y0_dbl_(y0),
        theta_(theta),
        theta_dbl_(value_of(theta)),
        // theta_fvars_(theta_dbl_.begin(), theta_dbl_.end()),
        x_(x),
        x_int_(x_int),
        N_(y0.size()),
        M_(theta.size()),
        size_(N_ + N_ * M_),
        msgs_(msgs) {}

  /**
   * Assign the derivative vector with the system derivatives at
   * the specified state and time.
   *
   * <p>The input state must be of size <code>size()</code>, and
   * the output produced will be of the same size.
   *
   * @param[in] z state of the coupled ode system.
   * @param[out] dz_dt populated with the derivatives of
   * the coupled system at the specified state and time.
   * @param[in]  t time.
   * @throw exception if the system function does not return the
   * same number of derivatives as the state vector size.
   *
   * y is the base ODE system state
   *
   */
  void operator()(const std::vector<double>& z, std::vector<double>& dz_dt,
                  double t) const {
    using std::vector;

    Eigen::Map<vector_d> dy_dt(dz_dt.data(), N_);
    matrix_d Jy_store(N_, N_);
    Eigen::Map<matrix_d> Jy(Jy_store.data(), N_, N_);
    //matrix_d Jtheta(N_, M_);
    

    //Eigen::Map<vector_d>(dz_dt.data(), N_) = dy_dt;

    // handle sensitivities wrt to theta
    Eigen::Map<const matrix_d> Stheta(z.data() + N_, N_, M_);
    Eigen::Map<matrix_d> dStheta_dt(dz_dt.data() + N_, N_, M_);

    //dStheta_dt = Jy * Stheta + Jtheta;
    jacobian_ode_states_parameters(f_, z, t, theta_dbl_, x_, x_int_, msgs_,
                                   dy_dt, Jy, dStheta_dt);
    dStheta_dt.noalias() += Jy * Stheta;
  }

  /**
   * Returns the size of the coupled system.
   *
   * @return size of the coupled system.
   */
  size_t size() const { return size_; }

  /**
   * Returns the initial state of the coupled system.  Because the
   * initial values are known, the initial state of the coupled
   * system is the same as the initial state of the base ODE
   * system.
   *
   * <p>This initial state returned is of size <code>size()</code>
   * where the first N (base ODE system size) parameters are the
   * initial conditions of the base ode system and the rest of the
   * initial condition elements are 0.
   *
   * @return the initial condition of the coupled system.
   */
  std::vector<double> initial_state() const {
    std::vector<double> state(size_, 0.0);
    for (size_t n = 0; n < N_; n++)
      state[n] = y0_dbl_[n];
    return state;
  }

  /**
   * Returns the base ODE system state corresponding to the
   * specified coupled system state.
   *
   * @param y coupled states after solving the ode
   */
  std::vector<std::vector<var>> decouple_states(
      const std::vector<std::vector<double>>& y) const {
    std::vector<var> temp_vars(N_);
    std::vector<double> temp_gradients(M_);
    std::vector<std::vector<var>> y_return(y.size());

    for (size_t i = 0; i < y.size(); i++) {
      // iterate over number of equations
      for (size_t j = 0; j < N_; j++) {
        // iterate over parameters for each equation
        for (size_t k = 0; k < M_; k++)
          temp_gradients[k] = y[i][y0_dbl_.size() + y0_dbl_.size() * k + j];

        temp_vars[j] = precomputed_gradients(y[i][j], theta_, temp_gradients);
      }
      y_return[i] = temp_vars;
    }
    return y_return;
  }
};

/**
 * The coupled ODE system for unknown initial values and known
 * parameters.
 *
 * <p>If the original ODE has states of size N, the
 * coupled system has N + N * N states. (derivatives of each
 * state with respect to each initial value)
 *
 * <p>The coupled system has N + N * N states, where N is the size of
 * the state vector in the base system.
 *
 * <p>The first N states correspond to the base system's N states:
 * \f$ \frac{d x_n}{dt} \f$
 *
 * <p>The next N states correspond to the sensitivities of the initial
 * conditions with respect to the to the first base system equation:
 * \f[
 *  \frac{d x_{N+n}}{dt}
 *     = \frac{d}{dt} \frac{\partial x_1}{\partial y0_n}
 * \f]
 *
 * <p>The next N states correspond to the sensitivities with respect
 * to the second base system equation, etc.
 *
 * @tparam F type of base ODE system functor
 */
template <typename F>
struct coupled_ode_system<F, var, double> {
  const F& f_;
  const std::vector<var>& y0_;
  const std::vector<double> y0_dbl_;
  const std::vector<double>& theta_dbl_;
  const std::vector<double>& x_;
  const std::vector<int>& x_int_;
  std::ostream* msgs_;
  const size_t N_;
  const size_t M_;
  const size_t size_;

  /**
   * Construct a coupled ODE system for an unknown initial state
   * and known parameters givne the specified base system functor,
   * base initial state, parameters, data, and an output stream
   * for messages.
   *
   * @param[in] f base ODE system functor.
   * @param[in] y0 initial state of the base ODE.
   * @param[in] theta system parameters.
   * @param[in] x real data.
   * @param[in] x_int integer data.
   * @param[in, out] msgs output stream for messages.
   */
  coupled_ode_system(const F& f, const std::vector<var>& y0,
                     const std::vector<double>& theta,
                     const std::vector<double>& x,
                     const std::vector<int>& x_int, std::ostream* msgs)
      : f_(f),
        y0_(y0),
        y0_dbl_(value_of(y0)),
        theta_dbl_(theta),
        x_(x),
        x_int_(x_int),
        msgs_(msgs),
        N_(y0.size()),
        M_(theta.size()),
        size_(N_ + N_ * N_) {}

  /**
   * Calculates the derivative of the coupled ode system
   * with respect to the state y at time t.
   *
   * @param[in] z the current state of the coupled, shifted ode
   * system. This is a a vector of double of length size().
   * @param[out] dz_dt a vector of length size() with the
   * derivatives of the coupled system evaluated with state y and
   * time t.
   * @param[in] t time.
   * @throw exception if the system functor does not return a
   * derivative vector of the same size as the state vector.
   *
   * y is the base ODE system state
   *
   */
  void operator()(const std::vector<double>& z, std::vector<double>& dz_dt,
                  double t) const {
    using std::vector;

    Eigen::Map<vector_d> dy_dt(dz_dt.data(), N_);

    matrix_d Jy_store(N_, N_);
    Eigen::Map<matrix_d> Jy(Jy_store.data(), N_, N_);

    jacobian_ode_states(f_, z, t, theta_dbl_, x_, x_int_, msgs_, dy_dt, Jy);

    // handle sensitivities wrt to theta
    Eigen::Map<const matrix_d> Sy(z.data() + N_, N_, N_);
    Eigen::Map<matrix_d> dSy_dt(dz_dt.data() + N_, N_, N_);

    dSy_dt.noalias() = Jy * Sy;
  }

  /**
   * Returns the size of the coupled system.
   *
   * @return size of the coupled system.
   */
  size_t size() const { return size_; }

  /**
   * Returns the initial state of the coupled system.
   *
   * <p>Because the starting state is unknown, the coupled system
   * incorporates the initial conditions as parameters.  The
   * initial conditions for the coupled part of the system are set
   * to zero along with the rest of the initial state, because the
   * value of the initial state has been moved into the
   * parameters.
   *
   * @return the initial condition of the coupled system.
   *   This is a vector of length size() where all elements
   *   are 0.
   */
  std::vector<double> initial_state() const {
    std::vector<double> initial(size_, 0.0);
    for (size_t i = 0; i < N_; i++)
      initial[i] = y0_dbl_[i];
    for (size_t i = 0; i < N_; i++)
      initial[N_ + i * N_ + i] = 1.0;
    return initial;
  }

  /**
   * Return the solutions to the basic ODE system, including
   * appropriate autodiff partial derivatives, given the specified
   * coupled system solution.
   *
   * @param y the vector of the coupled states after solving the ode
   */
  std::vector<std::vector<var>> decouple_states(
      const std::vector<std::vector<double>>& y) const {
    using std::vector;

    vector<var> temp_vars(N_);
    vector<double> temp_gradients(N_);
    vector<vector<var>> y_return(y.size());

    for (size_t i = 0; i < y.size(); i++) {
      // iterate over number of equations
      for (size_t j = 0; j < N_; j++) {
        // iterate over parameters for each equation
        for (size_t k = 0; k < N_; k++)
          temp_gradients[k] = y[i][y0_.size() + y0_.size() * k + j];

        temp_vars[j] = precomputed_gradients(y[i][j], y0_, temp_gradients);
      }
      y_return[i] = temp_vars;
    }

    return y_return;
  }
};

/**
 * The coupled ode system for unknown intial values and unknown
 * parameters.
 *
 * <p>The coupled system has N + N * (N + M) states, where N is
 * size of the base ODE state vector and M is the number of
 * parameters.
 *
 * <p>The first N states correspond to the base system's N states:
 *   \f$ \frac{d x_n}{dt} \f$
 *
 * <p>The next N+M states correspond to the sensitivities of the
 * initial conditions, then to the parameters with respect to the
 * to the first base system equation:
 *
 * \f[
 *   \frac{d x_{N + n}}{dt}
 *     = \frac{d}{dt} \frac{\partial x_1}{\partial y0_n}
 * \f]
 *
 * \f[
 *   \frac{d x_{N+N+m}}{dt}
 *     = \frac{d}{dt} \frac{\partial x_1}{\partial \theta_m}
 * \f]
 *
 * <p>The next N+M states correspond to the sensitivities with
 * respect to the second base system equation, etc.
 *
 * <p>If the original ode has a state vector of size N states and
 * a parameter vector of size M, the coupled system has N + N * (N
 * + M) states. (derivatives of each state with respect to each
 * initial value and each theta)
 *
 * @tparam F the functor for the base ode system
 */
template <typename F>
struct coupled_ode_system<F, var, var> {
  const F& f_;
  const std::vector<var>& y0_;
  const std::vector<double> y0_dbl_;
  const std::vector<var>& theta_;
  const std::vector<double> theta_dbl_;
  // mutable std::vector<fvar<double>> theta_fvars_;
  const std::vector<double>& x_;
  const std::vector<int>& x_int_;
  const size_t N_;
  const size_t M_;
  const size_t size_;
  std::ostream* msgs_;

  /**
   * Construct a coupled ODE system with unknown initial value and
   * known parameters, given the base ODE system functor, the
   * initial state of the base ODE, the parameters, data, and an
   * output stream to which to write messages.
   *
   * @param[in] f the base ode system functor.
   * @param[in] y0 the initial state of the base ode.
   * @param[in] theta parameters of the base ode.
   * @param[in] x real data.
   * @param[in] x_int integer data.
   * @param[in, out] msgs output stream to which to print messages.
   */
  coupled_ode_system(const F& f, const std::vector<var>& y0,
                     const std::vector<var>& theta,
                     const std::vector<double>& x,
                     const std::vector<int>& x_int, std::ostream* msgs)
      : f_(f),
        y0_(y0),
        y0_dbl_(value_of(y0)),
        theta_(theta),
        theta_dbl_(value_of(theta)),
        // theta_fvars_(theta_dbl_.begin(), theta_dbl_.end()),
        x_(x),
        x_int_(x_int),
        N_(y0.size()),
        M_(theta.size()),
        size_(N_ + N_ * (N_ + M_)),
        msgs_(msgs) {}

  /**
   * Populates the derivative vector with derivatives of the
   * coupled ODE system state with respect to time evaluated at the
   * specified state and specified time.
   *
   * @param[in]  z the current state of the coupled, shifted ode system,
   * of size <code>size()</code>.
   * @param[in, out] dz_dt populate with the derivatives of the
   * coupled system evaluated at the specified state and time.
   * @param[in] t time.
   * @throw exception if the base system does not return a
   * derivative vector of the same size as the state vector.
   *
   * y is the base ODE system state
   *
   */
  void operator()(const std::vector<double>& z, std::vector<double>& dz_dt,
                  double t) const {
    using std::vector;

    Eigen::Map<vector_d> dy_dt(dz_dt.data(), N_);

    Eigen::Map<const matrix_d> Sy(z.data() + N_, N_, N_);
    Eigen::Map<matrix_d> dSy_dt(dz_dt.data() + N_, N_, N_);

    Eigen::Map<const matrix_d> Stheta(z.data() + N_ + N_ * N_, N_, M_);
    Eigen::Map<matrix_d> dStheta_dt(dz_dt.data() + N_ + N_ * N_, N_, M_);

    /*
    jacobian_ode_states_parameters(f_, z, t, theta_dbl_, x_, x_int_, msgs_,
                                   dy_dt, Jy, Jtheta);

    // handle sensitivities wrt to y
    dSy_dt = Jy * Sy;

    // handle sensitivities wrt to theta
    dStheta_dt = Jy * Stheta + Jtheta;
    */

    // by reference value version
    // we write Jy into dSy_dt & Jtheta into dStheta_dt
    jacobian_ode_states_parameters(f_, z, t, theta_dbl_, x_, x_int_, msgs_,
                                   dy_dt, dSy_dt, dStheta_dt);

    dStheta_dt.noalias() += dSy_dt * Stheta;
    dSy_dt *= Sy;

    /* by value version with copying
    matrix_d Jy_store(N_, N_);
    Eigen::Map<matrix_d> Jy(Jy_store.data(), N_, N_);
    matrix_d Jtheta_store(N_, M_);
    Eigen::Map<matrix_d> Jtheta(Jtheta_store.data(), N_, M_);
    
    jacobian_ode_states_parameters(f_, z, t, theta_dbl_, x_, x_int_, msgs_,
                                   dy_dt, Jy, Jtheta);


    dStheta_dt.noalias() = Jy * Stheta + Jtheta;
    dSy_dt.noalias() = Jy * Sy;
    */

  }

  /**
   * Returns the size of the coupled system.
   *
   * @return size of the coupled system.
   */
  size_t size() const { return size_; }

  /**
   * Returns the initial state of the coupled system.
   *
   * Because the initial state is unknown, the coupled system
   * incorporates the initial condition offset from zero as
   * a parameter, and hence the return of this function is a
   * vector of zeros.
   *
   * @return the initial condition of the coupled system.  This is
   * a vector of length size() where all elements are 0.
   */
  std::vector<double> initial_state() const {
    std::vector<double> initial(size_, 0.0);
    for (size_t i = 0; i < N_; i++)
      initial[i] = y0_dbl_[i];
    for (size_t i = 0; i < N_; i++)
      initial[N_ + i * N_ + i] = 1.0;
    return initial;
  }

  /**
   * Return the basic ODE solutions given the specified coupled
   * system solutions, including the partials versus the
   * parameters encoded in the autodiff results.
   *
   * @param y the vector of the coupled states after solving the ode
   */
  std::vector<std::vector<var>> decouple_states(
      const std::vector<std::vector<double>>& y) const {
    using std::vector;

    vector<var> vars = y0_;
    vars.insert(vars.end(), theta_.begin(), theta_.end());

    vector<var> temp_vars(N_);
    vector<double> temp_gradients(N_ + M_);
    vector<vector<var>> y_return(y.size());

    for (size_t i = 0; i < y.size(); i++) {
      // iterate over number of equations
      for (size_t j = 0; j < N_; j++) {
        // iterate over parameters for each equation
        for (size_t k = 0; k < N_ + M_; k++)
          temp_gradients[k] = y[i][N_ + N_ * k + j];

        temp_vars[j] = precomputed_gradients(y[i][j], vars, temp_gradients);
      }
      y_return[i] = temp_vars;
    }

    return y_return;
  }
};

}  // namespace math
}  // namespace stan
#endif
