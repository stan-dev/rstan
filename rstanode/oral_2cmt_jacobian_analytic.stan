
functions {

  // slice out a sub-range from the vector also make one for int
  // arrays, etc.
  real[] slice_r(real[] rdata, int[] S, int s) {
    int l = sum(S[1:(s-1)]);
    int u = l + S[s];
    return(rdata[(l+1):u]);
  }
  int[] slice_i(int[] idata, int[] S, int s) {
    int l = sum(S[1:(s-1)]);
    int u = l + S[s];
    return(idata[(l+1):u]);
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

  // create an integer sequence
  int[] seq_int(int start, int end) {
    int N = end - start + 1;
    int seq[N];
    for(i in 1:N) seq[i] = i + start - 1;
    return(seq);
  }


  //ode-parameters: ka, ke, k12, k21
  //ode: da/dt = -ka * a
  //ode: dm/dt = ka * a - (k12 + ke) * m + k21 * p
  //ode: dp/dt = k12 * m - k21 * p
  //ode-pre-include:
  //ode-post-include:
  #include "oral_2cmt_jacobian_analytic_ode.stan"
  // note: the above stan file is automatically generated and must be
  // included. The filename is for foo.stan foo_ode.stan

  // this example does not use pre or post includes, but should
  // definitions needed to be included into the ODE RHS before the
  // dydt definition, then that must be put in separate files. Note
  // that there usually needs to be a pre.stan and a matching
  // pre.hpp. This is useful for forcing functions like dosing in PK
  // systems.

  real[,] integrate_serial(real[,] y0, real[] t0,
                           real[] ts, int[] Si_ts,
                           real[,] theta,
                           real[] x_r, int[] Si_x_r,
                           int[] x_i, int[] Si_x_i,
                           real rel_tol, real abs_tol, int max_steps,
                           int stiff) {
    int J = size(Si_ts)-1;

    real res[size(ts),3];
    for(j in 1:J) {
      if(stiff) {
        res[Si_ts[j] : Si_ts[j+1] - 1] = integrate_ode_bdf(oral_2cmt_jacobian_analytic_ode, y0[j], t0[j], ts[Si_ts[j] : Si_ts[j+1] - 1], theta[j], x_r[Si_x_r[j] : Si_x_r[j+1] - 1], x_i[Si_x_i[j] : Si_x_i[j+1] - 1],
                                                           rel_tol, abs_tol, max_steps);
      } else {
        res[Si_ts[j] : Si_ts[j+1] - 1] = integrate_ode_rk45(oral_2cmt_jacobian_analytic_ode, y0[j], t0[j], ts[Si_ts[j] : Si_ts[j+1] - 1], theta[j], x_r[Si_x_r[j] : Si_x_r[j+1] - 1], x_i[Si_x_i[j] : Si_x_i[j+1] - 1],
                                                            rel_tol, abs_tol, max_steps);
      }
    }
    return(res);
  }
}
data {
  int<lower=1> T;
  int<lower=1> J;
  real theta[4];
  real<lower=0> pbar[3 + 4];
  real theta_sd[4];
  real dose;
  int<lower=0,upper=1> parallel;
  int<lower=0,upper=1> observed;
  vector[ observed == 1 ? J*T : 0 ] yobs;
  real<lower=0> rel_tol;
  real<lower=0> abs_tol;
  int<lower=0> max_steps;
  int<lower=0,upper=2> solver;
}
transformed data {
  real state0[J,3];
  real t0[J];
  real ts[T*J];
  int S_ts[J];
  int Si_ts[J+1];
  real x_r[J];
  int S_x_r[J];
  int Si_x_r[J+1];
  int x_i[0];
  int S_x_i[J];
  int Si_x_i[J+1];
  vector[J] yobs_T;
  int stiff;

  stiff = solver == 1 ? 1 : 0;

  for(j in 1:J) {
    state0[j,1] = 0;
    state0[j,2] = 0;
    state0[j,3] = 0;

    x_r[j] = log(dose + j - 1.);

    t0[j] = 0;
    S_ts[j] = T;
    S_x_r[j] = 1;
    S_x_i[j] = 0;

    ts[(j-1) * T + 1 : j * T] = to_array_1d(seq_int(1, T));
  }

  Si_ts  = make_slice_index(S_ts);
  Si_x_r = make_slice_index(S_x_r);
  Si_x_i = make_slice_index(S_x_i);

  yobs_T = rep_vector(100., J);

  if(stiff) {
    print("Using stiff BDF integrator for serial integrator.");
  } else {
    print("Using non-stiff RK45 integrator for serial integrator.");
  }
  if(solver == 0) {
    print("Using non-stiff RK45 integrator for parallel integrator.");
  } else if(solver == 1) {
    print("Using stiff BDF integrator for parallel integrator.");
  } else if(solver == 2) {
    print("Using non-stiff Adams-Moulton integrator for parallel integrator.");
  }
  print("Relative tolerance: ", rel_tol);
  print("Absolute tolerance: ", abs_tol);
  print("Maximum # of steps: ", max_steps);

  if(parallel) {
    reject("Parallel integration not supported for now => go with map_rect.");
    print("Parallel ODE integration using OpenMP.");
  } else {
    print("Serial ODE integration.");
  }
}
parameters {
  real<lower=0> theta_v[4];
  real<lower=0> dose0_v[J];
  real<lower=0> sigma_y;
}
transformed parameters {
  //real yhat_par[sum(M),3];
  //real yhat_ser[sum(M),3];
  real yhat[sum(S_ts),3];
  real<lower=0> state0_v[J,3];
  real Theta_v[J,4];

  for(j in 1:J) {
    state0_v[j,1] = dose0_v[j];
    state0_v[j,2] = 0;
    state0_v[j,3] = 0;

    Theta_v[j] = theta_v;
  }

  if(parallel) {
    reject("Parallel integration not supported for now => go with map_rect.");
    //yhat = integrate_parallel(state0_v, t0, ts, Si_ts, Theta_v, x_r, Si_x_r, x_i, Si_x_i, rel_tol, abs_tol, max_steps, solver, pbar);
  } else {
    yhat = integrate_serial(state0_v, t0, ts, Si_ts, Theta_v, x_r, Si_x_r, x_i, Si_x_i, rel_tol, abs_tol, max_steps, stiff);
  }

}
model {
  vector[J] conc_T;

  for(j in 1:J)
    conc_T[j] = yhat[(j-1) * T + T,2];

  target += lognormal_lpdf(theta_v| log(theta), theta_sd);

  target += lognormal_lpdf(dose0_v| log(10.), 1);

  target += normal_lpdf(sigma_y| 0, 1);

  if(observed)
    target += lognormal_lpdf(yobs | log( to_vector(yhat[:,2]) + 1E-6 ), sigma_y);
  else
    target += normal_lpdf(yobs_T | conc_T, 5);
}
generated quantities {
  real yrepl[T*J];

  for(i in 1:(T*J))
    yrepl[i] = lognormal_rng(log(yhat[i,2] + 1E-6), sigma_y);
}
