#ifndef STAN_MATH_PRIM_ARR_FUNCTOR_integrate_1d_HPP
#define STAN_MATH_PRIM_ARR_FUNCTOR_integrate_1d_HPP

#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <functional>
#include <limits>
#include <ostream>
#include <vector>

namespace stan {
namespace math {
/**
 * Integrate a single variable function f from a to b to within a specified
 * tolerance. This function assumes a is less than b.
 *
 * The signature for f should be:
 *   double foo(double x, double xc)
 *
 * It should return the value of the function evaluated at x.
 *
 * Depending on whether or not a is finite or negative infinity and b is finite
 * or positive infinity, a different version of the 1d quadrature algorithm from
 * the Boost quadrature library is chosen.
 *
 * Integrals that cross zero are broken into two, and the separate integrals are
 * each integrated to the given tolerance.
 *
 * For integrals with finite limits, the xc argument is the distance to the
 * nearest boundary. So for a > 0, b > 0, it will be a - x for x closer to a,
 * and b - x for x closer to b. xc is computed in a way that avoids the
 * precision loss of computing a - x or b - x manually. For integrals that cross
 * zero, xc can take values a - x, -x, or b - x depending on which integration
 * limit it is nearest.
 *
 * If either limit is infinite, xc is set to NaN
 *
 * @tparam T Type of f
 * @param f the function to be integrated
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param tolerance target tolerance passed to Boost quadrature
 * @return numeric integral of function f
 */
template <typename F>
inline double integrate(const F& f, double a, double b, double tolerance) {
  double error1 = 0.0;
  double error2 = 0.0;
  double L1 = 0.0;
  double L2 = 0.0;
  bool used_two_integrals = false;
  size_t levels;
  double Q = 0.0;
  if (std::isinf(a) && std::isinf(b)) {
    auto f_wrap = [&](double x) {
      return f(x, std::numeric_limits<double>::quiet_NaN());
    };
    boost::math::quadrature::sinh_sinh<double> integrator;
    Q = integrator.integrate(f_wrap, tolerance, &error1, &L1, &levels);
  } else if (std::isinf(a)) {
    boost::math::quadrature::exp_sinh<double> integrator;
    /**
     * If the integral crosses zero, break it into two (advice from the Boost
     * implementation:
     * https://www.boost.org/doc/libs/1_66_0/libs/math/doc/html/math_toolkit/double_exponential/de_caveats.html)
     */
    if (b <= 0.0) {
      auto f_wrap = [&](double x) {
        return f(-(x + b), std::numeric_limits<double>::quiet_NaN());
      };
      Q = integrator.integrate(f_wrap, tolerance, &error1, &L1, &levels);
    } else {
      boost::math::quadrature::tanh_sinh<double> integrator_right;
      auto f_wrap = [&](double x) {
        return f(-x, std::numeric_limits<double>::quiet_NaN());
      };
      Q = integrator.integrate(f_wrap, tolerance, &error1, &L1, &levels)
          + integrator_right.integrate(f_wrap, -b, 0, tolerance, &error2, &L2,
                                       &levels);
      used_two_integrals = true;
    }
  } else if (std::isinf(b)) {
    boost::math::quadrature::exp_sinh<double> integrator;
    if (a >= 0.0) {
      auto f_wrap = [&](double x) {
        return f(x + a, std::numeric_limits<double>::quiet_NaN());
      };
      Q = integrator.integrate(f_wrap, tolerance, &error1, &L1, &levels);
    } else {
      boost::math::quadrature::tanh_sinh<double> integrator_right;
      auto f_wrap = [&](double x) {
        return f(x, std::numeric_limits<double>::quiet_NaN());
      };
      Q = integrator.integrate(f_wrap, tolerance, &error1, &L1, &levels)
          + integrator_right.integrate(f_wrap, a, 0, tolerance, &error2, &L2,
                                       &levels);
      used_two_integrals = true;
    }
  } else {
    auto f_wrap = [&](double x, double xc) { return f(x, xc); };
    boost::math::quadrature::tanh_sinh<double> integrator;
    if (a < 0.0 && b > 0.0) {
      Q = integrator.integrate(f_wrap, a, 0.0, tolerance, &error1, &L1, &levels)
          + integrator.integrate(f_wrap, 0.0, b, tolerance, &error2, &L2,
                                 &levels);
      used_two_integrals = true;
    } else {
      Q = integrator.integrate(f_wrap, a, b, tolerance, &error1, &L1, &levels);
    }
  }

  static const char* function = "integrate";
  if (used_two_integrals) {
    if (error1 > tolerance * L1) {
      domain_error(function, "error estimate of integral below zero", error1,
                   "",
                   " exceeds the given tolerance times integral below zero");
    }
    if (error2 > tolerance * L2) {
      domain_error(function, "error estimate of integral above zero", error2,
                   "",
                   " exceeds the given tolerance times integral above zero");
    }
  } else {
    if (error1 > tolerance * L1) {
      domain_error(function, "error estimate of integral", error1, "",
                   " exceeds the given tolerance times integral");
    }
  }
  return Q;
}

/**
 * Compute the integral of the single variable function f from a to b to within
 * a specified tolerance. a and b can be finite or infinite.
 *
 * The signature for f should be:
 *   double foo(double x, double xc, std::vector<double> theta,
 * std::vector<double> x_r, std::vector<int> x_i, std::ostream& msgs)
 *
 * It should return the value of the function evaluated at x. Any errors
 * should be printed to the msgs stream.
 *
 * Integrals that cross zero are broken into two, and the separate integrals are
 * each integrated to the given tolerance.
 *
 * For integrals with finite limits, the xc argument is the distance to the
 * nearest boundary. So for a > 0, b > 0, it will be a - x for x closer to a,
 * and b - x for x closer to b. xc is computed in a way that avoids the
 * precision loss of computing a - x or b - x manually. For integrals that cross
 * zero, xc can take values a - x, -x, or b - x depending on which integration
 * limit it is nearest.
 *
 * If either limit is infinite, xc is set to NaN
 *
 * Boost decides the integration is converged when subsequent estimates of the
 * integral are less than tolerance * abs(integral). This means the tolerance is
 * relative to the actual scale of the integral. Integrals that cross zero are
 * split into two. In this case, each integral is separately integrated to the
 * given tolerance.
 *
 * @tparam T Type of f
 * @param f the function to be integrated
 * @param a lower limit of integration
 * @param b upper limit of integration
 * @param theta additional parameters to be passed to f
 * @param x_r additional data to be passed to f
 * @param x_i additional integer data to be passed to f
 * @param[in, out] msgs the print stream for warning messages
 * @param tolerance integrator tolerance passed to Boost quadrature
 * @return numeric integral of function f
 */
template <typename F>
inline double integrate_1d(
    const F& f, const double a, const double b,
    const std::vector<double>& theta, const std::vector<double>& x_r,
    const std::vector<int>& x_i, std::ostream& msgs,
    const double tolerance
    = std::sqrt(std::numeric_limits<double>::epsilon())) {
  static const char* function = "integrate_1d";
  check_less_or_equal(function, "lower limit", a, b);

  if (a == b) {
    if (std::isinf(a))
      domain_error(function, "Integration endpoints are both", a, "", "");
    return 0.0;
  } else {
    return integrate(
        std::bind<double>(f, std::placeholders::_1, std::placeholders::_2,
                          theta, x_r, x_i, std::ref(msgs)),
        a, b, tolerance);
  }
}

}  // namespace math
}  // namespace stan

#endif
