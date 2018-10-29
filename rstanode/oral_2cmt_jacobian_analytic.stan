
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

  real[] oral_2cmt_ode(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    real dydt[3];
    real a = y[1];
    real m = y[2];
    real p = y[3];
    real ka = theta[1];
    real ke = theta[2];
    real k12 = theta[3];
    real k21 = theta[4];

    dydt[1] = -ka * a;
    dydt[2] = ka * a - (k12 + ke) * m + k21 * p;
    dydt[3] = k12 * m - k21 * p;

    return(dydt);
  }
  
  vector integrate_subject(vector theta, vector eta, real[] x_r, int[] x_i) {
    int T = x_i[1];
    int solver = x_i[2];
    int max_steps = x_i[3];
    int analytic_jacobian = x_i[4];
    real y0[3] = {10.0, 0, 0};
    real phi[4] = {theta[1], exp(log(theta[2]) + 0.2 * eta[1]), theta[3], theta[3]/theta[4]};
    int rel_tol_idx = 2*T+1;
    int abs_tol_idx = 2*T+2;
    real yhat[T,3];

    if(analytic_jacobian) {
      if(solver == 0) {
        yhat = integrate_ode_rk45(oral_2cmt_jacobian_analytic_ode, y0, 0.0, x_r[1:T], phi, x_r[1:1], x_i[1:1],
                                  x_r[rel_tol_idx], x_r[abs_tol_idx], max_steps);
      } else if(solver == 1) {
        yhat = integrate_ode_bdf(oral_2cmt_jacobian_analytic_ode, y0, 0.0, x_r[1:T], phi, x_r[1:1], x_i[1:1],
                                 x_r[rel_tol_idx], x_r[abs_tol_idx], max_steps);
      } else {
        yhat = integrate_ode_adams(oral_2cmt_jacobian_analytic_ode, y0, 0.0, x_r[1:T], phi, x_r[1:1], x_i[1:1],
                                   x_r[rel_tol_idx], x_r[abs_tol_idx], max_steps);
      }
    } else {
      if(solver == 0) {
        yhat = integrate_ode_rk45(oral_2cmt_ode, y0, 0.0, x_r[1:T], phi, x_r[1:1], x_i[1:1],
                                  x_r[rel_tol_idx], x_r[abs_tol_idx], max_steps);
      } else if(solver == 1) {
        yhat = integrate_ode_bdf(oral_2cmt_ode, y0, 0.0, x_r[1:T], phi, x_r[1:1], x_i[1:1],
                                 x_r[rel_tol_idx], x_r[abs_tol_idx], max_steps);
      } else {
        yhat = integrate_ode_adams(oral_2cmt_ode, y0, 0.0, x_r[1:T], phi, x_r[1:1], x_i[1:1],
                                   x_r[rel_tol_idx], x_r[abs_tol_idx], max_steps);
      }
    }

    return to_vector(yhat[:,2]);
  }

  vector subject_logLik(vector theta, vector eta, real[] x_r, int[] x_i) {
    int T = x_i[1];
    real sigma_y = theta[5];
    vector[T] log_yhat = log(integrate_subject(theta, eta, x_r, x_i) + 1E-6);

    return [ lognormal_lpdf(x_r[T+1:2*T] | log_yhat, sigma_y) ]';
  }
}
data {
  int<lower=0> T;
  int<lower=1> J;
  vector<lower=0>[4] prior_theta_mean;
  vector<lower=0>[4] prior_theta_sd;
  vector<lower=0>[J*T] yobs;
  real<lower=0> rel_tol;
  real<lower=0> abs_tol;
  int<lower=0> max_steps;
  int<lower=0,upper=2> solver;
  int<lower=0,upper=1> analytic_jacobian;
}
transformed data {
  int X_i[J,4];
  real X_r[J,2*T+2];

  for(j in 1:J) {
    X_i[j,1] = T;
    X_i[j,2] = solver;
    X_i[j,3] = max_steps;
    X_i[j,4] = analytic_jacobian;
    
    X_r[j,1:T] = to_array_1d(seq_int(1, T));
    X_r[j,(T+1):2*T] = to_array_1d(yobs[(j-1)*T+1:j*T]);
    X_r[j,2*T+1] = rel_tol;
    X_r[j,2*T+2] = abs_tol;
  }

  if(solver == 0) {
    print("Using non-stiff RK45 integrator.");
  } else if(solver == 1) {
    print("Using stiff BDF integrator.");
  } else if(solver == 2) {
    print("Using non-stiff Adams-Moulton integrator.");
  }

  if(analytic_jacobian) {
    print("Using analytic Jacobians for ODE integration.");
  } else {
    print("Using autodiff Jacobians for ODE integration.");
  }
  
  print("Relative tolerance: ", rel_tol);
  print("Absolute tolerance: ", abs_tol);
  print("Maximum # of steps: ", max_steps);
}
parameters {
  // ka, CL, Q, V2 (V1==1)
  vector<lower=0>[4] theta;
  vector[1] eta[J];
  real<lower=0> sigma_y;
}
transformed parameters {
}
model {
  target += lognormal_lpdf(theta| log(prior_theta_mean), prior_theta_sd);
  target += multi_normal_cholesky_lpdf(eta| [0]', [ [ 1 ] ]);
  target += normal_lpdf(sigma_y| 0, 0.5);
  target += sum(map_rect(subject_logLik, append_row(theta, sigma_y), eta, X_r, X_i));
}
generated quantities {
  vector[T*J] yrepl;
  vector[T*J] yhat = map_rect(integrate_subject, theta, eta, X_r, X_i);

  for(i in 1:(T*J))
    yrepl[i] = lognormal_rng(log(yhat[i] + 1E-6), sigma_y);
}
