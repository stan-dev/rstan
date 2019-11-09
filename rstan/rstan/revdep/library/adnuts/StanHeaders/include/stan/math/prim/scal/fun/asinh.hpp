#ifndef STAN_MATH_PRIM_SCAL_FUN_ASINH_HPP
#define STAN_MATH_PRIM_SCAL_FUN_ASINH_HPP

#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <stan/math/prim/scal/meta/likely.hpp>
#include <stan/math/prim/scal/fun/boost_policy.hpp>
#include <boost/math/special_functions/asinh.hpp>

namespace stan {
namespace math {

/**
 * Return the inverse hyperbolic sine of the specified value.
 * Returns infinity for infinity argument and -infinity for
 * -infinity argument.
 * Returns nan for nan argument.
 *
 * @param[in] x Argument.
 * @return Inverse hyperbolic sine of the argument.
 */
inline double asinh(double x) {
  if (unlikely(is_nan(x)))
    return x;
  else
    return boost::math::asinh(x, boost_policy_t());
}

/**
 * Integer version of asinh.
 *
 * @param[in] x Argument.
 * @return Inverse hyperbolic sine of the argument.
 */
inline double asinh(int x) { return asinh(static_cast<double>(x)); }

}  // namespace math
}  // namespace stan
#endif
