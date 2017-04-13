namespace stan {
  namespace math {
    template<>
    struct ode_system<ode_functor> {
      typedef ode_functor F;
      const F& f_;
      const std::vector<double> theta_;
      const std::vector<double>& x_r_;
      const std::vector<int>& x_i_;
      std::ostream* msgs_;

      ode_system(const F& f,
                 const std::vector<double> theta,
                 const std::vector<double>& x_r,
                 const std::vector<int>& x_i,
                 std::ostream* msgs)
	: f_(f),
	  theta_(theta),
	  x_r_(x_r),
	  x_i_(x_i),
	  msgs_(msgs)
      {
        //std::cout << "Yeah, this is the fast one!" << std::endl;
      }

      inline void operator()(const double t, const std::vector<double>& y,
                             std::vector<double>& dy_dt) const {
	using Eigen::VectorXd;
	using std::vector;
	using std::pow;
	using std::exp;

        Eigen::Map<VectorXd> dydt(&dy_dt[0], dy_dt.size());

#include ODE_HEADER(defs)
#include ODE_HEADER(rhs)

        // avoid call to functor which prevents a copy
        //dy_dt = f_(t, y, theta_, x_r_, x_i_, msgs_);
      }

      /*
      template <typename Derived1>
      inline void operator()(const double t, const std::vector<double>& y,
                             Eigen::MatrixBase<Derived1>& dydt) const {
	using Eigen::VectorXd;
	using std::vector;
	using std::pow;
	using std::exp;

#include ODE_HEADER(defs)
#include ODE_HEADER(rhs)

        // avoid call to functor which prevents a copy
        //dy_dt = f_(t, y, theta_, x_r_, x_i_, msgs_);
      }
      */

      template <typename Derived1, typename Derived2>
      inline void
      jacobian(const double t,
               const std::vector<double>& y,
               Eigen::MatrixBase<Derived1>& dydt,
               Eigen::MatrixBase<Derived2>& Jy
               ) const {
	using Eigen::VectorXd;
	using std::vector;
	using std::pow;
	using std::exp;

       	//const std::vector<double> f = f_(t, y, theta_, x_r_, x_i_, msgs_);
	//dydt = VectorXd::Map(&f[0], f.size());

#include ODE_HEADER(defs)
#include ODE_HEADER(rhs)
#include ODE_HEADER(jac_states)
      }

      template <typename Derived1, typename Derived2>
      inline void
      jacobian(const double t,
               const std::vector<double>& y,
               Eigen::MatrixBase<Derived1>& dydt,
               Eigen::MatrixBase<Derived2>& Jy,
               Eigen::MatrixBase<Derived2>& Jtheta
               ) const {
	using Eigen::VectorXd;
	using std::vector;
	using std::pow;
	using std::exp;

       	//const std::vector<double> f = f_(t, y, theta_, x_r_, x_i_, msgs_);
	//dydt = VectorXd::Map(&f[0], f.size());

#include ODE_HEADER(defs)
#include ODE_HEADER(rhs)
#include ODE_HEADER(jac_states)
#include ODE_HEADER(jac_params)
      }

    };

  }
}
