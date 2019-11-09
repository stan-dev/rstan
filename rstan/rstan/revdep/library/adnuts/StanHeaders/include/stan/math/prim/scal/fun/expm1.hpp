#ifndef STAN_MATH_PRIM_SCAL_FUN_EXPM1_HPP
#define STAN_MATH_PRIM_SCAL_FUN_EXPM1_HPP

#include <stan/math/prim/scal/fun/boost_policy.hpp>
#include <boost/math/special_functions/expm1.hpp>

namespace stan {
namespace math {

/**
 * Return the natural exponentiation of x minus one.
 * Returns infinity for infinity argument and -infinity for
 * -infinity argument.
 *
 * @param[in] x Argument.
 * @return Natural exponentiation of argument minus one.
 */
inline double expm1(double x) {
  return boost::math::expm1(x, boost_policy_t());
}

/**
 * Integer version of expm1.
 *
 * @param[in] x Argument.
 * @return Natural exponentiation of argument minus one.
 */
inline double expm1(int x) { return expm1(static_cast<double>(x)); }

}  // namespace math
}  // namespace stan
#endif
