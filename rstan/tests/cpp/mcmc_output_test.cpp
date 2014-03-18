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
    int argc = 0;
    char *argv[0];
    RInside R(argc, argv);
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

TEST_F(RStan, mcmc_writer) {
  std::stringstream sample_stream;
  std::stringstream diagnostic_stream;

  stan::io::mcmc_writer<stan_model> mcmc_writer(&stan_ss, &stan_ds);
  // Stan's usage:
  //   <nothing else>
  
  // then:
  //   writer.print_sample_names(s, sampler_ptr, model);
  //   writer.print_diagnostic_names(s, sampler_ptr, model);

  
  EXPECT_EQ("", stan_ss.str());
  EXPECT_EQ("", stan_ds.str());
}

