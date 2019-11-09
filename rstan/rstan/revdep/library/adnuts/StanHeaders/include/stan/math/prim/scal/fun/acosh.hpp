#ifndef STAN_MATH_PRIM_SCAL_FUN_ACOSH_HPP
#define STAN_MATH_PRIM_SCAL_FUN_ACOSH_HPP

#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <stan/math/prim/scal/meta/likely.hpp>
#include <stan/math/prim/scal/fun/boost_policy.hpp>
#include <boost/math/special_functions/acosh.hpp>

namespace stan {
namespace math {

/**
 * Return the inverse hyperbolic cosine of the specified value.
 * Returns nan for nan argument.
 *
 * @param[in] x Argument.
 * @return Inverse hyperbolic cosine of the argument.
 * @throw std::domain_error If argument is less than 1.
 */
inline double acosh(double x) {
  if (unlikely(is_nan(x)))
    return x;
  else
    return boost::math::acosh(x, boost_policy_t());
}

/**
 * Integer version of acosh.
 *
 * @param[in] x Argument.
 * @return Inverse hyperbolic cosine of the argument.
 * @throw std::domain_error If argument is less than 1.
 */
inline double acosh(int x) { return acosh(static_cast<double>(x)); }

}  // namespace math
}  // namespace stan
#endif
