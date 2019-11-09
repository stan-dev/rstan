#ifndef STAN_MATH_PRIM_SCAL_FUN_ERF_HPP
#define STAN_MATH_PRIM_SCAL_FUN_ERF_HPP

#include <stan/math/prim/scal/fun/boost_policy.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <limits>

namespace stan {
namespace math {

/**
 * Return the error function of the specified value.
 *
 * \f[
 * \mbox{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} dt
 * \f]
 *
 * @param[in] x Argument.
 * @return Error function of the argument.
 */
inline double erf(double x) {
  if (is_nan(x))
    return std::numeric_limits<double>::quiet_NaN();
  return boost::math::erf(x, boost_policy_t());
}

/**
 * Return the error function of the specified argument.  This
 * version is required to disambiguate <code>erf(int)</code>.
 *
 * @param[in] x Argument.
 * @return Error function of the argument.
 */
inline double erf(int x) { return erf(static_cast<double>(x)); }

}  // namespace math
}  // namespace stan
#endif
