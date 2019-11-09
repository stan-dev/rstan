#ifndef STAN_MATH_PRIM_SCAL_ERR_CHECK_POSITIVE_SIZE_HPP
#define STAN_MATH_PRIM_SCAL_ERR_CHECK_POSITIVE_SIZE_HPP

#include <stan/math/prim/scal/err/invalid_argument.hpp>
#include <sstream>
#include <string>

namespace stan {
namespace math {

/**
 * Check if <code>size</code> is positive.
 *
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param expr Expression for the dimension size (for error messages)
 * @param size Size value to check
 *
 * @throw <code>std::invalid_argument</code> if <code>size</code> is
 *   zero or negative.
 */
inline void check_positive_size(const char* function, const char* name,
                                const char* expr, int size) {
  if (size <= 0) {
    std::stringstream msg;
    msg << "; dimension size expression = " << expr;
    std::string msg_str(msg.str());
    invalid_argument(function, name, size, "must have a positive size, but is ",
                     msg_str.c_str());
  }
}

}  // namespace math
}  // namespace stan
#endif
