#ifndef STAN_MATH_PRIM_SCAL_FUN_TGAMMA_HPP
#define STAN_MATH_PRIM_SCAL_FUN_TGAMMA_HPP

#include <boost/math/special_functions/gamma.hpp>

namespace stan {
namespace math {

/**
 * Return the gamma function applied to the specified argument.
 *
 * @param x Argument.
 * @return The gamma function applied to argument.
 */
inline double tgamma(double x) { return boost::math::tgamma(x); }

}  // namespace math
}  // namespace stan
#endif
