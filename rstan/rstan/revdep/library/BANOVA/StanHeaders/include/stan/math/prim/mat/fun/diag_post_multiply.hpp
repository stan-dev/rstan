#ifndef STAN_MATH_PRIM_MAT_FUN_DIAG_POST_MULTIPLY_HPP
#define STAN_MATH_PRIM_MAT_FUN_DIAG_POST_MULTIPLY_HPP

#include <stan/math/prim/mat/err/check_vector.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <boost/math/tools/promotion.hpp>

namespace stan {
namespace math {

template <typename T1, typename T2, int R1, int C1, int R2, int C2>
Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type, R1, C1>
diag_post_multiply(const Eigen::Matrix<T1, R1, C1>& m1,
                   const Eigen::Matrix<T2, R2, C2>& m2) {
  check_vector("diag_post_multiply", "m2", m2);
  int m1_cols = m1.cols();
  check_size_match("diag_post_multiply", "m2.size()", m2.size(), "m1.cols()",
                   m1_cols);
  int m1_rows = m1.rows();
  Eigen::Matrix<typename boost::math::tools::promote_args<T1, T2>::type, R1, C1>
      result(m1_rows, m1_cols);

  for (int j = 0; j < m1_cols; ++j)
    for (int i = 0; i < m1_rows; ++i)
      result(i, j) = m2(j) * m1(i, j);
  return result;
}

}  // namespace math
}  // namespace stan
#endif
