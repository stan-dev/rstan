#ifndef STAN_MATH_GPU_ERR_CHECK_DIAGONAL_ZEROS_HPP
#define STAN_MATH_GPU_ERR_CHECK_DIAGONAL_ZEROS_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>

namespace stan {
namespace math {
/**
 * Check if the <code>matrix_gpu</code> has zeros on the diagonal
 *
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y <code>matrix_gpu</code> to test
 *
 * @throw <code>std::domain_error</code> if
 *    any diagonal element of the matrix is zero.
 */
inline void check_diagonal_zeros(const char* function, const char* name,
                                 const matrix_gpu& y) {
  if (y.size() == 0)
    return;

  cl::Kernel kernel_check_diagonal_zeros
      = opencl_context.get_kernel("is_zero_on_diagonal");
  cl::CommandQueue cmd_queue = opencl_context.queue();
  cl::Context ctx = opencl_context.context();
  try {
    int zero_on_diagonal_flag = 0;
    cl::Buffer buffer_flag(ctx, CL_MEM_READ_WRITE, sizeof(int));
    cmd_queue.enqueueWriteBuffer(buffer_flag, CL_TRUE, 0, sizeof(int),
                                 &zero_on_diagonal_flag);

    kernel_check_diagonal_zeros.setArg(0, y.buffer());
    kernel_check_diagonal_zeros.setArg(1, y.rows());
    kernel_check_diagonal_zeros.setArg(2, y.cols());
    kernel_check_diagonal_zeros.setArg(3, buffer_flag);

    cmd_queue.enqueueNDRangeKernel(kernel_check_diagonal_zeros, cl::NullRange,
                                   cl::NDRange(y.rows(), y.cols()),
                                   cl::NullRange);

    cmd_queue.enqueueReadBuffer(buffer_flag, CL_TRUE, 0, sizeof(int),
                                &zero_on_diagonal_flag);
    //  if zeros were found on the diagonal
    if (zero_on_diagonal_flag) {
      domain_error(function, name, "has zeros on the diagonal.", "");
    }
  } catch (const cl::Error& e) {
    check_opencl_error("diag_zeros_check", e);
  }
}

}  // namespace math
}  // namespace stan
#endif
#endif
