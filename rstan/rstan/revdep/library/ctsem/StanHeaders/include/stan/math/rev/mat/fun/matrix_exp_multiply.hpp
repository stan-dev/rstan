#ifndef STAN_MATH_REV_MAT_FUN_MATRIX_EXP_MULTIPLY_HPP
#define STAN_MATH_REV_MAT_FUN_MATRIX_EXP_MULTIPLY_HPP

#include <stan/math/rev/mat.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/matrix_exp_action_handler.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/mat/fun/to_var.hpp>
#include <stan/math/prim/mat/fun/matrix_exp.hpp>
#include <stan/math/prim/mat/err/check_multiplicable.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Calculate adjoint of matrix exponential action
 * exp(At)*B when A is data and B is var, used for chaining action.
 *
 * @param Ad double array pointer to the data of matrix A
 * @param n dim of square matrix A
 * @param adjexpAB MatrixXd adjoint of exp(At)*B
 * @param t double time data
 * @return MatrixXd The adjoint of B.
 */
inline Eigen::MatrixXd exp_action_chain_dv(double* Ad, const int& n,
                                           const Eigen::MatrixXd& adjexpAB,
                                           const double t) {
  using Eigen::Map;
  matrix_exp_action_handler handle;
  return handle.action(Map<Eigen::MatrixXd>(Ad, n, n).transpose(), adjexpAB, t);
}

/**
 * Calculate adjoint of matrix exponential action
 * exp(At)*B when A is var and B is data, used for chaining action.
 *
 * @param Ad double array pointer to the data of matrix A
 * @param Bd double array pointer to the data of matrix B
 * @param n dim(nb. of rows) of square matrix A
 * @param m nb. of cols of matrix B
 * @param adjexpAB MatrixXd adjoint of exp(At)*B
 * @param t double time data
 * @return MatrixXd The adjoint of A.
 */
inline Eigen::MatrixXd exp_action_chain_vd(double* Ad, double* Bd, const int& n,
                                           const int& m,
                                           const Eigen::MatrixXd& adjexpAB,
                                           const double t) {
  using Eigen::Map;
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  Eigen::MatrixXd adjexpA = Eigen::MatrixXd::Zero(n, n);
  Eigen::MatrixXd adjA = Eigen::MatrixXd::Zero(n, n);

  // TODO(yizhang): a better way such as complex step approximation
  try {
    start_nested();

    adjexpA = adjexpAB * Map<MatrixXd>(Bd, n, m).transpose();
    Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> Av(n, n);
    for (int i = 0; i < Av.size(); ++i) {
      Av(i) = to_var(Ad[i]);
    }
    std::vector<stan::math::var> Avec(Av.data(), Av.data() + Av.size());
    Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> expA
        = matrix_exp(Av);
    std::vector<double> g;
    for (size_type i = 0; i < expA.size(); ++i) {
      stan::math::set_zero_all_adjoints_nested();
      expA.coeffRef(i).grad(Avec, g);
      for (size_type j = 0; j < adjA.size(); ++j) {
        adjA(j) += adjexpA(i) * g[j];
      }
    }
  } catch (const std::exception& e) {
    recover_memory_nested();
    throw;
  }
  recover_memory_nested();
  return adjA;
}

/**
 * This is a subclass of the vari class for matrix
 * exponential action exp(At) * B where A is a double
 * NxN matrix and B is a NxCb matrix.
 *
 * The class stores the structure of each matrix,
 * the double values of A and B, and pointers to
 * the varis for A and B if A or B is a var. It
 * also instantiates and stores pointers to
 * varis for all elements of A * B.
 *
 * @tparam Ta Scalar type of matrix A
 * @tparam N rows and cols of A
 * @tparam Tb Scalar type for matrix B
 * @tparam Cb cols for matrix B
 */
template <typename Ta, int N, typename Tb, int Cb>
class matrix_exp_action_vari : public vari {
 public:
  int n_;
  int B_cols_;
  int A_size_;
  int B_size_;
  double t_;
  double* Ad_;
  double* Bd_;
  vari** variRefA_;
  vari** variRefB_;
  vari** variRefexpAB_;

  /**
   * Constructor: vari child-class of matrix_exp_action_vari.
   * @param A statically-sized matirx
   * @param B statically-sized matirx
   * @param t double scalar time.
   */
  matrix_exp_action_vari(const Eigen::Matrix<Ta, N, N>& A,
                         const Eigen::Matrix<Tb, N, Cb>& B, const double& t)
      : vari(0.0),
        n_(A.rows()),
        B_cols_(B.cols()),
        A_size_(A.size()),
        B_size_(B.size()),
        t_(t),
        Ad_(ChainableStack::instance().memalloc_.alloc_array<double>(A_size_)),
        Bd_(ChainableStack::instance().memalloc_.alloc_array<double>(B_size_)),
        variRefA_(
            ChainableStack::instance().memalloc_.alloc_array<vari*>(A_size_)),
        variRefB_(
            ChainableStack::instance().memalloc_.alloc_array<vari*>(B_size_)),
        variRefexpAB_(
            ChainableStack::instance().memalloc_.alloc_array<vari*>(B_size_)) {
    using Eigen::Map;
    using Eigen::MatrixXd;
    for (size_type i = 0; i < A.size(); ++i) {
      variRefA_[i] = A.coeffRef(i).vi_;
      Ad_[i] = A.coeffRef(i).val();
    }
    for (size_type i = 0; i < B.size(); ++i) {
      variRefB_[i] = B.coeffRef(i).vi_;
      Bd_[i] = B.coeffRef(i).val();
    }
    matrix_exp_action_handler handle;
    MatrixXd expAB = handle.action(Map<MatrixXd>(Ad_, n_, n_),
                                   Map<MatrixXd>(Bd_, n_, B_cols_), t_);
    for (size_type i = 0; i < expAB.size(); ++i)
      variRefexpAB_[i] = new vari(expAB.coeffRef(i), false);
  }

  /**
   * Chain command for the adjoint of matrix exp action.
   */
  virtual void chain() {
    using Eigen::Map;
    using Eigen::MatrixXd;
    MatrixXd adjexpAB(n_, B_cols_);

    for (size_type i = 0; i < adjexpAB.size(); ++i)
      adjexpAB(i) = variRefexpAB_[i]->adj_;

    MatrixXd adjA = exp_action_chain_vd(Ad_, Bd_, n_, B_cols_, adjexpAB, t_);
    MatrixXd adjB = exp_action_chain_dv(Ad_, n_, adjexpAB, t_);

    for (size_type i = 0; i < A_size_; ++i) {
      variRefA_[i]->adj_ += adjA(i);
    }
    for (size_type i = 0; i < B_size_; ++i) {
      variRefB_[i]->adj_ += adjB(i);
    }
  }
};

/**
 * This is a subclass of the vari class for matrix
 * exponential action exp(At) * B where A is an N by N
 * matrix of double, B is N by Cb, and t is data(time)
 *
 * The class stores the structure of each matrix,
 * the double values of A and B, and pointers to
 * the varis for A and B if B is a var. It
 * also instantiates and stores pointers to
 * varis for all elements of exp(At) * B.
 *
 * @tparam N Rows and cols for matrix A, also rows for B
 * @tparam Tb Scalar type for matrix B
 * @tparam Cb Columns for matrix B
 */
template <int N, typename Tb, int Cb>
class matrix_exp_action_vari<double, N, Tb, Cb> : public vari {
 public:
  int n_;
  int B_cols_;
  int A_size_;
  int B_size_;
  double t_;
  double* Ad_;
  double* Bd_;
  vari** variRefB_;
  vari** variRefexpAB_;

  /**
   * Constructor for matrix_exp_action_vari.
   *
   * All memory allocated in
   * ChainableStack's stack_alloc arena.
   *
   * It is critical for the efficiency of this object
   * that the constructor create new varis that aren't
   * popped onto the var_stack_, but rather are
   * popped onto the var_nochain_stack_. This is
   * controlled to the second argument to
   * vari's constructor.
   *
   * @param A matrix
   * @param B matrix
   * @param t double
   */
  matrix_exp_action_vari(const Eigen::Matrix<double, N, N>& A,
                         const Eigen::Matrix<Tb, N, Cb>& B, const double& t)
      : vari(0.0),
        n_(A.rows()),
        B_cols_(B.cols()),
        A_size_(A.size()),
        B_size_(B.size()),
        t_(t),
        Ad_(ChainableStack::instance().memalloc_.alloc_array<double>(A_size_)),
        Bd_(ChainableStack::instance().memalloc_.alloc_array<double>(B_size_)),
        variRefB_(
            ChainableStack::instance().memalloc_.alloc_array<vari*>(B_size_)),
        variRefexpAB_(
            ChainableStack::instance().memalloc_.alloc_array<vari*>(B_size_)) {
    using Eigen::Map;
    using Eigen::MatrixXd;
    for (size_type i = 0; i < A.size(); ++i)
      Ad_[i] = A.coeffRef(i);
    for (size_type i = 0; i < B.size(); ++i) {
      variRefB_[i] = B.coeffRef(i).vi_;
      Bd_[i] = B.coeffRef(i).val();
    }
    matrix_exp_action_handler handle;
    MatrixXd expAB = handle.action(Map<MatrixXd>(Ad_, n_, n_),
                                   Map<MatrixXd>(Bd_, n_, B_cols_), t_);
    for (size_type i = 0; i < expAB.size(); ++i)
      variRefexpAB_[i] = new vari(expAB.coeffRef(i), false);
  }

  /**
   * Chain for matrix_exp_action_vari.
   */
  virtual void chain() {
    using Eigen::Map;
    using Eigen::MatrixXd;
    MatrixXd adjexpAB(n_, B_cols_);
    // MatrixXd adjB(n_, B_cols_);

    for (size_type i = 0; i < adjexpAB.size(); ++i)
      adjexpAB(i) = variRefexpAB_[i]->adj_;

    MatrixXd adjB = exp_action_chain_dv(Ad_, n_, adjexpAB, t_);

    for (size_type i = 0; i < B_size_; ++i) {
      variRefB_[i]->adj_ += adjB(i);
    }
  }
};

/**
 * This is a subclass of the vari class for matrix
 * exponential action exp(At) * B where A is an N by N
 * matrix of var, B is N by Cb matrix of double, and t is data(time)
 *
 * The class stores the structure of each matrix,
 * the double values of A and B, and pointers to
 * the varis for A and B if B is a var. It
 * also instantiates and stores pointers to
 * varis for all elements of exp(At) * B.
 *
 * @tparam Ta Scalar type for matrix A
 * @tparam N Rows and cols for matrix A, also rows for B
 * @tparam Cb Columns for matrix B
 */
template <typename Ta, int N, int Cb>
class matrix_exp_action_vari<Ta, N, double, Cb> : public vari {
 public:
  int n_;
  int B_cols_;
  int A_size_;
  int B_size_;
  double t_;
  double* Ad_;
  double* Bd_;
  vari** variRefA_;
  vari** variRefexpAB_;

  /**
   * Constructor for matrix_exp_action_vari.
   *
   * All memory allocated in
   * ChainableStack's stack_alloc arena.
   *
   * It is critical for the efficiency of this object
   * that the constructor create new varis that aren't
   * popped onto the var_stack_, but rather are
   * popped onto the var_nochain_stack_. This is
   * controlled to the second argument to
   * vari's constructor.
   *
   * @param A matrix
   * @param B matrix
   * @param t double
   */
  matrix_exp_action_vari(const Eigen::Matrix<Ta, N, N>& A,
                         const Eigen::Matrix<double, N, Cb>& B, const double& t)
      : vari(0.0),
        n_(A.rows()),
        B_cols_(B.cols()),
        A_size_(A.size()),
        B_size_(B.size()),
        t_(t),
        Ad_(ChainableStack::instance().memalloc_.alloc_array<double>(A_size_)),
        Bd_(ChainableStack::instance().memalloc_.alloc_array<double>(B_size_)),
        variRefA_(
            ChainableStack::instance().memalloc_.alloc_array<vari*>(A_size_)),
        variRefexpAB_(
            ChainableStack::instance().memalloc_.alloc_array<vari*>(B_size_)) {
    using Eigen::Map;
    using Eigen::MatrixXd;
    for (size_type i = 0; i < A.size(); ++i) {
      variRefA_[i] = A.coeffRef(i).vi_;
      Ad_[i] = A.coeffRef(i).val();
    }
    for (size_type i = 0; i < B.size(); ++i) {
      Bd_[i] = B.coeffRef(i);
    }
    matrix_exp_action_handler handle;
    MatrixXd expAB = handle.action(Map<MatrixXd>(Ad_, n_, n_),
                                   Map<MatrixXd>(Bd_, n_, B_cols_), t_);
    for (size_type i = 0; i < expAB.size(); ++i)
      variRefexpAB_[i] = new vari(expAB.coeffRef(i), false);
  }

  virtual void chain() {
    using Eigen::Map;
    using Eigen::MatrixXd;
    MatrixXd adjexpAB(n_, B_cols_);

    for (size_type i = 0; i < adjexpAB.size(); ++i)
      adjexpAB(i) = variRefexpAB_[i]->adj_;

    MatrixXd adjA = exp_action_chain_vd(Ad_, Bd_, n_, B_cols_, adjexpAB, t_);

    for (size_type i = 0; i < A_size_; ++i) {
      variRefA_[i]->adj_ += adjA(i);
    }
  }
};

/**
 * Return product of exp(At) and B, where A is a NxN matrix,
 * B is a NxCb matrix, and t is a double
 * @tparam Ta scalar type matrix A
 * @tparam N Rows and cols matrix A, also rows of matrix B
 * @tparam Tb scalar type matrix B
 * @tparam Cb Columns matrix B
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @param[in] t double
 * @return exponential of At multiplies B
 */
template <typename Ta, int N, typename Tb, int Cb>
inline typename boost::enable_if_c<boost::is_same<Ta, var>::value
                                       || boost::is_same<Tb, var>::value,
                                   Eigen::Matrix<var, N, Cb> >::type
matrix_exp_action(const Eigen::Matrix<Ta, N, N>& A,
                  const Eigen::Matrix<Tb, N, Cb>& B, const double& t = 1.0) {
  matrix_exp_action_vari<Ta, N, Tb, Cb>* baseVari
      = new matrix_exp_action_vari<Ta, N, Tb, Cb>(A, B, t);
  Eigen::Matrix<var, N, Cb> expAB_v(A.rows(), B.cols());
  for (size_type i = 0; i < expAB_v.size(); ++i) {
    expAB_v.coeffRef(i).vi_ = baseVari->variRefexpAB_[i];
  }
  return expAB_v;
}

/**
 * Wrapper of matrix_exp_action function for a more literal name
 * @tparam Ta scalar type matrix A
 * @tparam Tb scalar type matrix B
 * @tparam Cb Columns matrix B
 * @param[in] A Matrix
 * @param[in] B Matrix
 * @return exponential of A multiplies B
 */
template <typename Ta, typename Tb, int Cb>
inline Eigen::Matrix<typename stan::return_type<Ta, Tb>::type, -1, Cb>
matrix_exp_multiply(const Eigen::Matrix<Ta, -1, -1>& A,
                    const Eigen::Matrix<Tb, -1, Cb>& B) {
  check_nonzero_size("matrix_exp_multiply", "input matrix", A);
  check_nonzero_size("matrix_exp_multiply", "input matrix", B);
  check_multiplicable("matrix_exp_multiply", "A", A, "B", B);
  check_square("matrix_exp_multiply", "input matrix", A);
  return matrix_exp_action(A, B);
}

}  // namespace math
}  // namespace stan

#endif
