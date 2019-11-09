R"(
#define A(i, j)  A[j * rows + i]
#define B(i, j)  B[j * rows + i]
__kernel void copy(
      __global double *A,
      __global double *B,
      unsigned int rows,
      unsigned int cols) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < rows && j < cols) {
    B(i, j) = A(i, j);
  }
};)"
