R"(
#ifndef A
#define A(i, j)  A[j * rows + i]
#endif
__kernel void is_nan(
      __global double *A,
      int rows,
      int cols,
      __global int *flag) {
  const int i = get_global_id(0);
  const int j = get_global_id(1);
  if (i < rows && j < cols) {
    if (isnan(A(i, j))) {
      flag[0] = 1;
    }
  }
};)"
