R"(
#ifndef A
#define A(i, j)  A[j * rows + i]
#endif
__kernel void is_symmetric(
      __global double *A,
      int rows,
      int cols,
      __global int *flag,
      double tolerance) {
  const int i = get_global_id(0);
  const int j = get_global_id(1);
  if (i < rows && j < cols) {
    double diff = fabs(A(i, j) - A(j, i));
    if (diff > tolerance) {
      flag[0] = 0;
    }
  }
};)"
