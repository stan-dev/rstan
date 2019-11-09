PARAMETER_SECTION

  vector prob(1,nobs)    // per capita mort prob

PROCEDURE_SECTION

  dvariable fpen=0.0;         // penalty variable
  // power-Ricker
  prob = ccc*pow(elem_prod(TBL/dddd,exp(1-TBL/dddd)),ggggggg);
  // penalties: constrain 0.001 <= prob <= 0.999
  prob = posfun(prob,0.001,fpen);
  f += 1000*fpen;
  prob = 1-posfun(1-prob,0.001,fpen);
  f += 1000*fpen;
  // binomial negative log-likelihood
  f -= sum( log_comb(nexposed,Kill)+
            elem_prod(Kill,log(prob))+
            elem_prod(nexposed-Kill,log(1-prob)));
