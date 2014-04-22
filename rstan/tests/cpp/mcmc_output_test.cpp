#include <gtest/gtest.h>
#include <RInside.h>
#include <rstan/stan_fit.hpp>  // for find_index
#include <rstan/mcmc_output.hpp>
#include <stan/io/mcmc_writer.hpp>
#include <sstream>
#include "blocker.cpp"

typedef boost::ecuyer1988 rng_t; // (2**50 = 1T samples, 1000 chains)
typedef stan::mcmc::adapt_dense_e_nuts<stan_model, rng_t> sampler;

//using namespace rstan;

class mock_sampler : public stan::mcmc::base_mcmc {
public:
  mock_sampler()
    : base_mcmc(&std::cout, &std::cout), n_transition_called(0) { }

  stan::mcmc::sample transition(stan::mcmc::sample& init_sample) {
    n_transition_called++;
    return init_sample;
  }

  int n_transition_called;
};

class RStan : public ::testing::Test {
public:
  static void SetUpTestCase() { 
    RInside R;
  }

  static void TearDownTestCase() { }
  
  void SetUp() { 

    std::fstream data_stream(std::string("tests/cpp/blocker.data.R").c_str());
    stan::io::dump data_context(data_stream);
    data_stream.close();
    
    model_ptr = new stan_model(data_context, &std::cout);

    cont_params = Eigen::VectorXd::Zero(model_ptr->num_params_r());

    sampler_ptr = new mock_sampler;
  }
  
  void TearDown() { 
    delete model_ptr;
    delete sampler_ptr;
  }
  


  std::stringstream ss, ds, stan_ss, stan_ds;
  stan_model* model_ptr;
  Eigen::VectorXd cont_params;
  mock_sampler* sampler_ptr;
};

TEST_F(RStan, mcmc_output) {
  stan::mcmc::sample s(cont_params, 0, 0);
  int ctrl_sampling_iter_save = 1000;

  rstan::mcmc_output<stan_model> mcmc_output(&ss, &ds);
  mcmc_output.inspect_internal("after construction");

  std::vector<std::string> iter_param_names;
  std::vector<std::string> sampler_param_names;
  // RStan's usage:
  
  mcmc_output.set_output_names(s, sampler_ptr, *model_ptr, 
                               iter_param_names, sampler_param_names);
  mcmc_output.inspect_internal("after set_output_names");
  std::cout << "iter_param_names: " << std::endl;
  for (int i = 0; i < iter_param_names.size(); i++) {
    std::cout << "  " << iter_param_names[i] << std::endl;
  }
  std::cout << "sampler_param_names: " << std::endl;
  for (int i = 0; i < sampler_param_names.size(); i++) {
    std::cout << "  " << sampler_param_names[i] << std::endl;
  }
  
  std::vector<Rcpp::NumericVector> sampler_params;
  mcmc_output.init_sampler_params(sampler_params, ctrl_sampling_iter_save);
  mcmc_output.inspect_internal("after init_sampler_params");
  std::cout << "sampler_params[" << sampler_params.size() << "]:" << std::endl;
  for (int i = 0; i < sampler_params.size(); i++)
    std::cout << "  " << sampler_params[i] << std::endl;

  std::vector<Rcpp::NumericVector> iter_params;
  mcmc_output.init_iter_params(iter_params, ctrl_sampling_iter_save);
  mcmc_output.inspect_internal("after init_iter_params");
  std::cout << "iter_params[" << iter_params.size() << "]:" << std::endl;
  for (int i = 0; i < iter_params.size(); i++)
    std::cout << "  " << iter_params[i] << std::endl;
  
  std::cout << "iter_params[0].size(): " << iter_params[0].size() << std::endl;

  ASSERT_EQ(2, iter_params.size());
  for (int i = 0; i < 2; i++)
    ASSERT_EQ(ctrl_sampling_iter_save, iter_params[i].size());
  
  EXPECT_EQ("", ss.str());
  EXPECT_EQ("", ds.str());

  ss.str("");
  ds.str("");
  mcmc_output.print_sample_names();
  EXPECT_EQ("lp__,accept_stat__,d,sigmasq_delta,mu.1,mu.2,mu.3,mu.4,mu.5,mu.6,mu.7,mu.8,mu.9,mu.10,mu.11,mu.12,mu.13,mu.14,mu.15,mu.16,mu.17,mu.18,mu.19,mu.20,mu.21,mu.22,delta.1,delta.2,delta.3,delta.4,delta.5,delta.6,delta.7,delta.8,delta.9,delta.10,delta.11,delta.12,delta.13,delta.14,delta.15,delta.16,delta.17,delta.18,delta.19,delta.20,delta.21,delta.22,delta_new,sigma_delta\n", ss.str());
  EXPECT_EQ("", ds.str());
  
  ss.str("");
  ds.str("");
  mcmc_output.output_diagnostic_names(s, sampler_ptr, *model_ptr);
  EXPECT_EQ("", ss.str());
  EXPECT_EQ("lp__,accept_stat__\n", ds.str());
}

TEST_F(RStan, mcmc_output_inside_run_markov_chain) {
  // setup
  double lp = 123;
  ASSERT_EQ(47, model_ptr->num_params_r());
  Eigen::VectorXd x(47);
  for (int i = 0; i < x.size(); i++) {
    x(i) = i + 1;
  }
  stan::mcmc::sample s(x, lp, 0);
  
  int ctrl_sampling_iter_save = 1000;
  std::vector<std::string> iter_param_names;
  std::vector<std::string> sampler_param_names;
  std::vector<Rcpp::NumericVector> sampler_params;
  std::vector<Rcpp::NumericVector> iter_params;

  rstan::mcmc_output<stan_model> mcmc_output(&ss, &ds);
  mcmc_output.set_output_names(s, sampler_ptr, *model_ptr, 
                               iter_param_names, sampler_param_names);
  mcmc_output.init_sampler_params(sampler_params, ctrl_sampling_iter_save);
  mcmc_output.init_iter_params(iter_params, ctrl_sampling_iter_save);


  ASSERT_EQ(2, iter_params.size());
  for (int i = 0; i < 2; i++)
    ASSERT_EQ(ctrl_sampling_iter_save, iter_params[i].size());
  ASSERT_EQ("", ss.str());
  ASSERT_EQ("", ds.str());

  ASSERT_EQ(0, sampler_params.size());
  ASSERT_EQ(2, iter_params.size());
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < ctrl_sampling_iter_save; j++) {
      EXPECT_FLOAT_EQ(0.0, iter_params[i][j]);
    }
  }

  SUCCEED() << "the mcmc_output object is created as it is passed into run_markov_chain()";
  //----------------------------------------
  // iteration 1.

  ss.str("");
  ds.str("");

  //mcmc_output.inspect_internal("");

  rng_t base_rng;
  std::vector<Rcpp::NumericVector> chains;
  for (int i = 0; i < 2; i++)
    chains.push_back(Rcpp::NumericVector(ctrl_sampling_iter_save));
  bool is_warmup = true;
  std::vector<double> sum_pars(47);
  double sum_lp = 0;
  std::vector<size_t> qoi_idx;
  qoi_idx.push_back(0);
  qoi_idx.push_back(1);
  int iter_save_i = 0;
  std::stringstream rcout;

  mcmc_output.output_sample_params(base_rng, s, sampler_ptr, *model_ptr,
                                   chains, is_warmup,
                                   sampler_params, iter_params,
                                   sum_pars, sum_lp, qoi_idx,
                                   iter_save_i, &rcout);
  EXPECT_FLOAT_EQ(0.0, sum_lp);
  ASSERT_EQ(47, sum_pars.size());
  for (int i = 0; i < 47; i++) {
    EXPECT_FLOAT_EQ(0.0, sum_pars[i]);
  }
  ASSERT_EQ(0, sampler_params.size());
  ASSERT_EQ(2, iter_params.size());
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < ctrl_sampling_iter_save; j++) {
      if (i == 0 && j == 0) {
        EXPECT_FLOAT_EQ(lp, iter_params[i][j]) 
          << "i=" << i << ", j=" << j;
      }
      else {
        EXPECT_FLOAT_EQ(0.0, iter_params[i][j]) 
          << "i=" << i << ", j=" << j;
      }
    }
  }
  EXPECT_EQ("", rcout.str());
  EXPECT_EQ("123,0,1,7.38906,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,2.71828\n", 
            ss.str());
  EXPECT_EQ("", ds.str());

  ss.str("");
  ds.str("");
  mcmc_output.output_diagnostic_params(s, sampler_ptr);
  EXPECT_EQ("", ss.str());
  EXPECT_EQ("123,0\n", ds.str());

  //------------------------------------------------------------
  // iteration 2

  is_warmup = false;
  iter_save_i = 1;
  rcout.str("");
  ss.str("");
  ds.str("");
  mcmc_output.output_sample_params(base_rng, s, sampler_ptr, *model_ptr,
                                   chains, is_warmup,
                                   sampler_params, iter_params,
                                   sum_pars, sum_lp, qoi_idx,
                                   iter_save_i, &rcout);
  EXPECT_FLOAT_EQ(1.0 * lp, sum_lp);
  ASSERT_EQ(47, sum_pars.size());
  for (int i = 0; i < 47; i++) {
    if (i == 1) {
      EXPECT_FLOAT_EQ(exp(x(i)), sum_pars[i]);
    } else {
      EXPECT_FLOAT_EQ(x(i), sum_pars[i]) << "index: " << i << std::endl;
    }
  }
  ASSERT_EQ(0, sampler_params.size());
  ASSERT_EQ(2, iter_params.size());
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < ctrl_sampling_iter_save; j++) {
      if (i == 0 && j <= iter_save_i) {
        EXPECT_FLOAT_EQ(lp, iter_params[i][j]) 
          << "i=" << i << ", j=" << j;
      }
      else {
        EXPECT_FLOAT_EQ(0.0, iter_params[i][j]) 
          << "i=" << i << ", j=" << j;
      }
    }
  }
  EXPECT_EQ("", rcout.str());
  EXPECT_EQ("123,0,1,7.38906,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,2.71828\n", 
            ss.str());
  EXPECT_EQ("", ds.str());

  ss.str("");
  ds.str("");
  mcmc_output.output_diagnostic_params(s, sampler_ptr);
  EXPECT_EQ("", ss.str());
  EXPECT_EQ("123,0\n", ds.str());
  
  //------------------------------------------------------------
  // iteration 3

  is_warmup = false;
  iter_save_i = 2;
  rcout.str("");
  ss.str("");
  ds.str("");
  mcmc_output.output_sample_params(base_rng, s, sampler_ptr, *model_ptr,
                                   chains, is_warmup,
                                   sampler_params, iter_params,
                                   sum_pars, sum_lp, qoi_idx,
                                   iter_save_i, &rcout);

  EXPECT_FLOAT_EQ(2.0 * lp, sum_lp);
  ASSERT_EQ(47, sum_pars.size());
  for (int i = 0; i < 47; i++) {
    if (i == 1) {
      EXPECT_FLOAT_EQ(2 * exp(x(i)), sum_pars[i]);
    } else {
      EXPECT_FLOAT_EQ(2 * x(i), sum_pars[i]) << "index: " << i << std::endl;
    }
  }
  ASSERT_EQ(0, sampler_params.size());
  ASSERT_EQ(2, iter_params.size());
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < ctrl_sampling_iter_save; j++) {
      if (i == 0 && j <= iter_save_i) {
        EXPECT_FLOAT_EQ(lp, iter_params[i][j]) 
          << "i=" << i << ", j=" << j;
      }
      else {
        EXPECT_FLOAT_EQ(0.0, iter_params[i][j]) 
          << "i=" << i << ", j=" << j;
      }
    }
  }
  EXPECT_EQ("", rcout.str());
  EXPECT_EQ("123,0,1,7.38906,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,2.71828\n", 
            ss.str());
  EXPECT_EQ("", ds.str());


  ss.str("");
  ds.str("");
  mcmc_output.output_diagnostic_params(s, sampler_ptr);
  EXPECT_EQ("", ss.str());
  EXPECT_EQ("123,0\n", ds.str());

}


TEST_F(RStan, mcmc_writer) {
  double lp = 123;
  ASSERT_EQ(47, model_ptr->num_params_r());
  Eigen::VectorXd x(47);
  for (int i = 0; i < x.size(); i++) {
    x(i) = i + 1;
  }
  stan::mcmc::sample s(x, lp, 0);
  
  stan::common::recorder::csv stan_sample_recorder(&stan_ss, "# ");
  stan::common::recorder::csv stan_diagnostic_recorder(&stan_ds, "# ");
  stan::common::recorder::messages stan_message_recorder(&stan_ss, "# ");
  stan::io::mcmc_writer<stan_model,
                        stan::common::recorder::csv,
                        stan::common::recorder::csv,
                        stan::common::recorder::messages> 
    mcmc_writer(stan_sample_recorder, 
                stan_diagnostic_recorder, 
                stan_message_recorder,
                &std::cout);
  // Stan's usage:
  //   <nothing else>
  
  // then:
  //   writer.write_sample_names(s, sampler_ptr, model);
  mcmc_writer.write_sample_names(s, sampler_ptr, *model_ptr);
  EXPECT_EQ("lp__,accept_stat__,d,sigmasq_delta,mu.1,mu.2,mu.3,mu.4,mu.5,mu.6,mu.7,mu.8,mu.9,mu.10,mu.11,mu.12,mu.13,mu.14,mu.15,mu.16,mu.17,mu.18,mu.19,mu.20,mu.21,mu.22,delta.1,delta.2,delta.3,delta.4,delta.5,delta.6,delta.7,delta.8,delta.9,delta.10,delta.11,delta.12,delta.13,delta.14,delta.15,delta.16,delta.17,delta.18,delta.19,delta.20,delta.21,delta.22,delta_new,sigma_delta\n", 
            stan_ss.str());
  EXPECT_EQ("", stan_ds.str());

  
  stan_ss.str("");
  stan_ds.str("");
  mcmc_writer.write_diagnostic_names(s, sampler_ptr, *model_ptr);

  EXPECT_EQ("", stan_ss.str());
  EXPECT_EQ("lp__,accept_stat__\n", stan_ds.str());
}

TEST_F(RStan, mcmc_writer_inside_run_markov_chain) {
  double lp = 123;
  ASSERT_EQ(47, model_ptr->num_params_r());
  Eigen::VectorXd x(47);
  for (int i = 0; i < x.size(); i++) {
    x(i) = i + 1;
  }
  stan::mcmc::sample s(x, lp, 0);
  stan::common::recorder::csv stan_sample_recorder(&stan_ss, "# ");
  stan::common::recorder::csv stan_diagnostic_recorder(&stan_ds, "# ");
  stan::common::recorder::messages stan_message_recorder(&stan_ss, "# ");
  stan::io::mcmc_writer<stan_model,
                        stan::common::recorder::csv,
                        stan::common::recorder::csv,
                        stan::common::recorder::messages> 
    mcmc_writer(stan_sample_recorder, 
                stan_diagnostic_recorder, 
                stan_message_recorder,
                &std::cout);
  mcmc_writer.write_sample_names(s, sampler_ptr, *model_ptr);
  mcmc_writer.write_diagnostic_names(s, sampler_ptr, *model_ptr);

  rng_t base_rng;
  
  //------------------------------------------------------------
  // iteration 1
  stan_ss.str("");
  stan_ds.str("");
  mcmc_writer.write_sample_params(base_rng, s, *sampler_ptr, *model_ptr);
  EXPECT_EQ("123,0,1,7.38906,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,2.71828\n", stan_ss.str());
  EXPECT_EQ("", stan_ds.str());
  
  stan_ss.str("");
  stan_ds.str("");
  mcmc_writer.write_diagnostic_params(s, sampler_ptr);
  EXPECT_EQ("", stan_ss.str());
  EXPECT_EQ("123,0\n", stan_ds.str());

  //------------------------------------------------------------
  // iteration 2
  stan_ss.str("");
  stan_ds.str("");
  mcmc_writer.write_sample_params(base_rng, s, *sampler_ptr, *model_ptr);
  EXPECT_EQ("123,0,1,7.38906,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,2.71828\n", stan_ss.str());
  EXPECT_EQ("", stan_ds.str());
  
  stan_ss.str("");
  stan_ds.str("");
  mcmc_writer.write_diagnostic_params(s, sampler_ptr);
  EXPECT_EQ("", stan_ss.str());
  EXPECT_EQ("123,0\n", stan_ds.str());


  //------------------------------------------------------------
  // iteration 3
  stan_ss.str("");
  stan_ds.str("");
  mcmc_writer.write_sample_params(base_rng, s, *sampler_ptr, *model_ptr);
  EXPECT_EQ("123,0,1,7.38906,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,2.71828\n", 
            stan_ss.str());
  EXPECT_EQ("", stan_ds.str());
  
  stan_ss.str("");
  stan_ds.str("");
  mcmc_writer.write_diagnostic_params(s, sampler_ptr);
  EXPECT_EQ("", stan_ss.str());
  EXPECT_EQ("123,0\n", stan_ds.str());
}

