#ifndef STAN_MATH_REV_MAT_FUNCTOR_ODEINT_ODE_SYSTEM_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_ODEINT_ODE_SYSTEM_HPP

#include <stan/math/prim/arr/meta/get.hpp>
#include <stan/math/prim/arr/meta/length.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
//#include <stan/math/prim/mat/functor/coupled_ode_system.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
//#include <stan/math/rev/mat/functor/ode_system.hpp>
#include <stan/math/rev/core.hpp>
#include <ostream>
#include <stdexcept>
#include <vector>

namespace stan {
  namespace math {

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
    class odeint_ode_system <F, double, stan::math::var> {
      const std::vector<double> y0_dbl_;
      const size_t N_;
      const size_t M_;
      const size_t S_;
      const size_t size_;
      const ode_system<F> ode_system_;

    public:
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
      odeint_ode_system(const F& f,
                         const std::vector<double>& y0,
                         const std::vector<stan::math::var>& theta,
                         const std::vector<double>& x,
                         const std::vector<int>& x_int,
                         std::ostream* msgs)
        : y0_dbl_(y0),
          N_(y0.size()),
          M_(theta.size()),
          S_(M_),
          size_(N_ + N_ * M_),
          ode_system_(f, stan::math::value_of(theta), x, x_int, msgs) {}

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
      void operator()(const std::vector<double>& z,
                      std::vector<double>& dz_dt,
                      double t) const {
        using Eigen::Map;
        using Eigen::MatrixXd;
        using Eigen::VectorXd;

        const std::vector<double> y_base(z.begin(), z.begin() + N_);
        //Map<const VectorXd> y_base(&z[0], N_);
        Map<VectorXd> dz_dt_eig(&dz_dt[0], N_);
        Map<MatrixXd> dZ_dt_sens(&dz_dt[0] + N_, N_, S_);
        Map<const MatrixXd> Z_sens(&z[0] + N_, N_, S_);

        MatrixXd Jy(N_, N_);
        MatrixXd Jtheta(N_, M_);

        ode_system_.jacobian(t, y_base, dz_dt_eig, Jy, Jtheta);

        dZ_dt_sens = Jy * Z_sens + Jtheta;
       }

      /**
       * Returns the size of the coupled system.
       *
       * @return size of the coupled system.
       */
      size_t size() const {
        return size_;
      }

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
    class odeint_ode_system <F, stan::math::var, double> {
      const std::vector<double> y0_dbl_;
      const size_t N_;
      const size_t M_;
      const size_t S_;
      const size_t size_;
      const ode_system<F> ode_system_;

    public:
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
      odeint_ode_system(const F& f,
                         const std::vector<stan::math::var>& y0,
                         const std::vector<double>& theta,
                         const std::vector<double>& x,
                         const std::vector<int>& x_int,
                         std::ostream* msgs)
        : y0_dbl_(stan::math::value_of(y0)),
          N_(y0.size()),
          M_(theta.size()),
          S_(N_),
          size_(N_ + N_ * N_),
          ode_system_(f, theta, x, x_int, msgs) { }

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
      void operator()(const std::vector<double>& z,
                      std::vector<double>& dz_dt,
                      double t) const {
        using Eigen::Map;
        using Eigen::MatrixXd;
        using Eigen::VectorXd;

        const std::vector<double> y_base(z.begin(), z.begin() + N_);
        //Map<const VectorXd> y_base(&z[0], N_);
        Map<VectorXd> dz_dt_eig(&dz_dt[0], N_);
        Map<MatrixXd> dZ_dt_sens(&dz_dt[0] + N_, N_, S_);
        Map<const MatrixXd> Z_sens(&z[0] + N_, N_, S_);

        MatrixXd Jy(N_, N_);

        ode_system_.jacobian(t, y_base, dz_dt_eig, Jy);

        dZ_dt_sens = Jy * Z_sens;
      }

      /**
       * Returns the size of the coupled system.
       *
       * @return size of the coupled system.
       */
      size_t size() const {
        return size_;
      }

      /**
       * Returns the initial state of the coupled system.
       *
       * <p>For unkown initial values, the sensitivities are
       * initialized to the identity matrix at t=0.
       *
       * @return the initial condition of the coupled system.
       */
      std::vector<double> initial_state() const {
        std::vector<double> initial(size_, 0.0);
        for (size_t i = 0; i < N_; i++)
          initial[i] = y0_dbl_[i];
        for (size_t i = 0; i < N_; i++)
          initial[N_ + i * N_ + i] = 1.0;
        return initial;
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
    class odeint_ode_system <F, stan::math::var, stan::math::var> {
      const std::vector<double> y0_dbl_;
      const size_t N_;
      const size_t M_;
      const size_t S_;
      const size_t size_;
      const ode_system<F> ode_system_;

    public:
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
      odeint_ode_system(const F& f,
                         const std::vector<stan::math::var>& y0,
                         const std::vector<stan::math::var>& theta,
                         const std::vector<double>& x,
                         const std::vector<int>& x_int,
                         std::ostream* msgs)
        : y0_dbl_(stan::math::value_of(y0)),
          N_(y0.size()),
          M_(theta.size()),
          S_(N_ + M_),
          size_(N_ + N_ * (N_ + M_)),
          ode_system_(f, stan::math::value_of(theta), x, x_int, msgs) { }

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
      void operator()(const std::vector<double>& z,
                      std::vector<double>& dz_dt,
                      double t) const {
        using Eigen::Map;
        using Eigen::MatrixXd;
        using Eigen::VectorXd;

        const std::vector<double> y_base(z.begin(), z.begin() + N_);
        //Map<const VectorXd> y_base(&z[0], N_);
        Map<VectorXd> dz_dt_eig(&dz_dt[0], N_);
        Map<MatrixXd> dZ_dt_sens(&dz_dt[0] + N_, N_, S_);
        Map<const MatrixXd> Z_sens(&z[0] + N_, N_, S_);
        // write Jtheta directly into correct position of dZ_dt
        //Map<MatrixXd> dZ_dt_sens_param(&dz_dt[0] + N_ + N_ * N_, N_, M_);

        MatrixXd Jy(N_, N_);
        MatrixXd Jtheta(N_, M_);

        ode_system_.jacobian(t, y_base, dz_dt_eig, Jy, Jtheta);
        //dZ_dt_sens_param = Jtheta;

        dZ_dt_sens.rightCols(M_) = Jtheta;
        dZ_dt_sens.leftCols(N_).setZero();
        dZ_dt_sens += Jy * Z_sens;
      }

      /**
       * Returns the size of the coupled system.
       *
       * @return size of the coupled system.
       */
      size_t size() const {
        return size_;
      }

      /**
       * Returns the initial state of the coupled system.
       *
       * Because the initial state is unknown, its initial is set to
       * the identity matrix.
       *
       * @return the initial condition of the coupled system.
       */
      std::vector<double> initial_state() const {
        std::vector<double> initial(size_, 0.0);
        for (size_t i = 0; i < N_; i++)
          initial[i] = y0_dbl_[i];
        for (size_t i = 0; i < N_; i++)
          initial[N_ + i * N_ + i] = 1.0;
        return initial;
      }
    };
  }  // math
}  // stan

#endif
