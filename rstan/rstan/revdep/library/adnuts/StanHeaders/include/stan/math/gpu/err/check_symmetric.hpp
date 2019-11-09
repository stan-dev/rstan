#ifndef STAN_MATH_GPU_ERR_CHECK_SYMMETRIC_HPP
#define STAN_MATH_GPU_ERR_CHECK_SYMMETRIC_HPP
#ifdef STAN_OPENCL
#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>

namespace stan {
namespace math {
/**
 * Check if the <code>matrix_gpu</code> is symmetric
 *
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y <code>matrix_gpu</code> to test
 *
 * @throw <code>std::domain_error</code> if
 *    the matrix is not symmetric.
 */
inline void check_symmetric(const char* function, const char* name,
                            const matrix_gpu& y) {
  if (y.size() == 0)
    return;
  check_square(function, name, y);
  cl::Kernel kernel_check_symmetric = opencl_context.get_kernel("is_symmetric");
  cl::CommandQueue cmd_queue = opencl_context.queue();
  cl::Context& ctx = opencl_context.context();
  try {
    int symmetric_flag = 1;
    cl::Buffer buffer_symmetric_flag(ctx, CL_MEM_READ_WRITE, sizeof(int));
    cmd_queue.enqueueWriteBuffer(buffer_symmetric_flag, CL_TRUE, 0, sizeof(int),
                                 &symmetric_flag);

    kernel_check_symmetric.setArg(0, y.buffer());
    kernel_check_symmetric.setArg(1, y.rows());
    kernel_check_symmetric.setArg(2, y.cols());
    kernel_check_symmetric.setArg(3, buffer_symmetric_flag);
    kernel_check_symmetric.setArg(4, math::CONSTRAINT_TOLERANCE);

    cmd_queue.enqueueNDRangeKernel(kernel_check_symmetric, cl::NullRange,
                                   cl::NDRange(y.rows(), y.cols()),
                                   cl::NullRange);

    cmd_queue.enqueueReadBuffer(buffer_symmetric_flag, CL_TRUE, 0, sizeof(int),
                                &symmetric_flag);
    //  if the matrix is not symmetric
    if (!symmetric_flag) {
      domain_error(function, name, "is not symmetric", "");
    }
  } catch (const cl::Error& e) {
    check_opencl_error("symmetric_check", e);
  }
}

}  // namespace math
}  // namespace stan
#endif
#endif
