PARAMETER_SECTION

  vector herdvec(1,nobs)
  vector eta(1,nobs)
  vector mu(1,nobs)

PROCEDURE_SECTION

  herdvec = sigma_herd*(Zherd*u_herd);
  eta = X*beta;                                                // form linear predictor 
  eta += herdvec;                                              // augment with random effects
  mu = pow(1.0+exp(-eta),-1.0);                                // logistic transform
  f -= sum(elem_prod(incidence,log(mu))+elem_prod(size-incidence,log(1.0-mu))); // binomial log-likelihood (unnormalized)
  
  f+=0.5*norm2(u_herd);  // log-prior (standard normal)

