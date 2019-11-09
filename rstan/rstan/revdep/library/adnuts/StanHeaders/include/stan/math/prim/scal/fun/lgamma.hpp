#ifndef STAN_MATH_PRIM_SCAL_FUN_LGAMMA_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LGAMMA_HPP

#include <stan/math/prim/scal/fun/boost_policy.hpp>
#include <boost/math/special_functions/gamma.hpp>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the gamma function applied to
 * the specified argument.
 *
   \f[
   \mbox{lgamma}(x) =
   \begin{cases}
     \textrm{error} & \mbox{if } x\in \{\dots, -3, -2, -1, 0\}\\
     \ln\Gamma(x) & \mbox{if } x\not\in \{\dots, -3, -2, -1, 0\}\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{lgamma}(x)}{\partial x} =
   \begin{cases}
     \textrm{error} & \mbox{if } x\in \{\dots, -3, -2, -1, 0\}\\
     \Psi(x) & \mbox{if } x\not\in \{\dots, -3, -2, -1, 0\}\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
\f]
*
* @param x argument
* @return natural logarithm of the gamma function applied to
* argument
*/
inline double lgamma(double x) {
  return boost::math::lgamma(x, boost_policy_t());
}

/**
 * Return the natural logarithm of the gamma function applied
 * to the specified argument.
 *
 * @param x argument
 * @return natural logarithm of the gamma function applied to
 * argument
 */
inline double lgamma(int x) { return boost::math::lgamma(x, boost_policy_t()); }

}  // namespace math
}  // namespace stan
#endif
