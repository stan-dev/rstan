#ifndef STAN_MATH_REV_MAT_FUN_SCALE_MATRIX_EXP_MULTIPLY_HPP
#define STAN_MATH_REV_MAT_FUN_SCALE_MATRIX_EXP_MULTIPLY_HPP

#include <stan/math/prim/mat.hpp>
#include <stan/math/rev/mat/fun/matrix_exp_multiply.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Wrapper of matrix_exp_action function for a more literal name.
 * @tparam Ta scalar type matrix A
 * @tparam Tb scalar type matrix B
 * @tparam Cb Columns matrix B
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @param[in] t double
 * @return exponential of At multiplies B
 */
template <typename Ta, typename Tb, int Cb>
inline Eigen::Matrix<typename stan::return_type<Ta, Tb>::type, -1, Cb>
scale_matrix_exp_multiply(const double& t, const Eigen::Matrix<Ta, -1, -1>& A,
                          const Eigen::Matrix<Tb, -1, Cb>& B) {
  check_nonzero_size("scale_matrix_exp_multiply", "input matrix", A);
  check_nonzero_size("scale_matrix_exp_multiply", "input matrix", B);
  check_multiplicable("scale_matrix_exp_multiply", "A", A, "B", B);
  check_square("scale_matrix_exp_multiply", "input matrix", A);
  return matrix_exp_action(A, B, t);
}

}  // namespace math
}  // namespace stan

#endif
