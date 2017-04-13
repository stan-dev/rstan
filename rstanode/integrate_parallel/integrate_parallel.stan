  // slice out a sub-range from the vector also make one for int
  // arrays, etc.
  real[] slice_r(real[] data, int[] S, int s) {
    int l = sum(S[1:(s-1)]);
    int u = l + S[s];
    return(data[(l+1):u]);
  }
  int[] slice_i(int[] data, int[] S, int s) {
    int l = sum(S[1:(s-1)]);
    int u = l + S[s];
    return(data[(l+1):u]);
  }

 // turn a slicing variable for a ragged array
  // S = {5, 6, 11}
  // into
  // Si = {0, 5, 5 + 6, 5+6+11} + 1
  // such that we can index the ragged array A as
  // A[Si[u] : Si[u+1]-1]
  // for the uth unit
  int[] make_slice_index(int[] S) {
    int Si[size(S)+1];
    int cv = 1;
    Si[1] = cv;
    for(i in 1:size(S)) {
      cv = cv + S[i];
      Si[i+1] = cv;
    }
    return(Si);
  }

  // externally defined function which executes in parallel the ODE
  // solver
  real[,] integrate_parallel(real[,] y0, real[] t0,
                             real[] ts, int[] Si_ts,
                             real[,] theta,
                             real[] x_r, int[] Si_x_r,
                             int[] x_i, int[] Si_x_i,
                             real rel_tol, real abs_tol, int max_steps,
                             int solver, // 0=RK45 non-stiff; =1 =>
                                          // bdf stiff; =2 => Adams
                                          // non-stiff
                             real[] pbar);

