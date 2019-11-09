DATA_SECTION
	init_number nt			//Number of timepoints	
	init_vector dvar(2,nt)		//Change in variance
	init_number max_var_g		//Upper bound
PARAMETER_SECTION
	init_bounded_number rho(0.00,1.00)
	init_bounded_number var_g(0.001,max_var_g)
	likeprof_number prof_rho
	number dt
	number rss
	number predict
	number resid
	vector t(2,nt)
	objective_function_value f
	
PRELIMINARY_CALCS_SECTION
 dt=1/(nt-1);
 for(int i=2; i<=nt; i++)
 {
	t[i]=i*dt-dt;
 }
	
PROCEDURE_SECTION
	rss=0.0;
	prof_rho=rho;
	
	for(int i=2; i<=nt; i++)
	{		
		predict=pow(rho,2)*pow(dt,2)*var_g+
                        (1-pow(rho,2))*dt*var_g+
                        2*pow(rho,2)*dt*var_g*(t[i]-dt);
		resid=dvar[i]-predict;
		rss+=pow(resid,2); 
	}	
	f=nt/2*log(rss/nt);


