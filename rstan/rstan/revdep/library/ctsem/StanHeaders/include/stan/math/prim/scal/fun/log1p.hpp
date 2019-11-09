#ifndef STAN_MATH_PRIM_SCAL_FUN_LOG1P_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LOG1P_HPP

#include <stan/math/prim/scal/fun/boost_policy.hpp>
#include <boost/math/special_functions/log1p.hpp>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of one plus the specified value.
 *
 * \f[
 * \mbox{log1p}(x) = \log(1 + x)
 * \f]
 *
 * This version is more stable for arguments near zero than
 * the direct definition.  If <code>log1p(x)</code> is defined to
 * be negative infinity.
 *
 * @param[in] x Argument.
 * @return Natural log of one plus the argument.
 * @throw std::domain_error If argument is less than -1.
 */
inline double log1p(double x) {
  return boost::math::log1p(x, boost_policy_t());
}

/**
 * Return the natural logarithm of one plus the specified
 * argument.  This version is required to disambiguate
 * <code>log1p(int)</code>.
 *
 * @param[in] x Argument.
 * @return Natural logarithm of one plus the argument.
 * @throw std::domain_error If argument is less than -1.
 */
inline double log1p(int x) { return log1p(static_cast<double>(x)); }

}  // namespace math
}  // namespace stan
#endif
