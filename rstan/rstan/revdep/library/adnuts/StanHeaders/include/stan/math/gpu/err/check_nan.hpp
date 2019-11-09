#ifndef STAN_MATH_GPU_ERR_CHECK_NAN_HPP
#define STAN_MATH_GPU_ERR_CHECK_NAN_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>

namespace stan {
namespace math {
/**
 * Check if the <code>matrix_gpu</code> has NaN values
 *
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y <code>matrix_gpu</code> to test
 *
 * @throw <code>std::domain_error</code> if
 *    any element of the matrix is <code>NaN</code>.
 */
inline void check_nan(const char* function, const char* name,
                      const matrix_gpu& y) {
  if (y.size() == 0)
    return;

  cl::Kernel kernel_check_nan = opencl_context.get_kernel("is_nan");
  cl::CommandQueue cmd_queue = opencl_context.queue();
  cl::Context& ctx = opencl_context.context();
  try {
    int nan_flag = 0;
    cl::Buffer buffer_nan_flag(ctx, CL_MEM_READ_WRITE, sizeof(int));
    cmd_queue.enqueueWriteBuffer(buffer_nan_flag, CL_TRUE, 0, sizeof(int),
                                 &nan_flag);

    kernel_check_nan.setArg(0, y.buffer());
    kernel_check_nan.setArg(1, y.rows());
    kernel_check_nan.setArg(2, y.cols());
    kernel_check_nan.setArg(3, buffer_nan_flag);

    cmd_queue.enqueueNDRangeKernel(kernel_check_nan, cl::NullRange,
                                   cl::NDRange(y.rows(), y.cols()),
                                   cl::NullRange);

    cmd_queue.enqueueReadBuffer(buffer_nan_flag, CL_TRUE, 0, sizeof(int),
                                &nan_flag);
    //  if NaN values were found in the matrix
    if (nan_flag) {
      domain_error(function, name, "has NaN values", "");
    }
  } catch (const cl::Error& e) {
    check_opencl_error("nan_check", e);
  }
}

}  // namespace math
}  // namespace stan
#endif
#endif
