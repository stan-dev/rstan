#ifndef STAN_MATH_PRIM_SCAL_FUN_ERFC_HPP
#define STAN_MATH_PRIM_SCAL_FUN_ERFC_HPP

#include <stan/math/prim/scal/fun/boost_policy.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <limits>

namespace stan {
namespace math {

/**
 * Return the complementary error function of the specified value.
 *
 * \f[
 * \mbox{erfc}(x) = 1 - \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} dt
 * \f]
 *
 * @param[in] x Argument.
 * @return Complementary error function of the argument.
 */
inline double erfc(double x) {
  if (is_nan(x))
    return std::numeric_limits<double>::quiet_NaN();
  return boost::math::erfc(x, boost_policy_t());
}

/**
 * Return the error function of the specified argument.  This
 * version is required to disambiguate <code>erf(int)</code>.
 *
 * @param[in] x Argument.
 * @return Complementary error function value of the argument.
 */
inline double erfc(int x) { return erfc(static_cast<double>(x)); }

}  // namespace math
}  // namespace stan
#endif
