R"(
#ifndef A
#define A(i, j)  A[j * rows + i]
#endif
__kernel void is_zero_on_diagonal(
        __global double *A,
        int rows,
        int cols,
        __global int *flag) {
  const int i = get_global_id(0);
  if (i < rows && i < cols) {
    if (A(i, i) == 0) {
      flag[0] = 1;
    }
  }
};)"
