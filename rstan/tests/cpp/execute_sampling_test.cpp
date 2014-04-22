#include <gtest/gtest.h>
#include <RInside.h>
#include <rstan/stan_fit.hpp>  
#include "blocker.cpp"
#include <sstream>
#include <fstream>
#include <boost/algorithm/string.hpp>

typedef boost::ecuyer1988 rng_t; // (2**50 = 1T samples, 1000 chains)
typedef stan::mcmc::adapt_dense_e_nuts<stan_model, rng_t> sampler_t;

class RStan : public ::testing::Test {
public:
  RStan() 
    : in_(),
      data_stream_(std::string("tests/cpp/blocker.data.R").c_str()),
      data_context_(data_stream_),
      model_(data_context_, &std::cout),
      sample_stream(),
      diagnostic_stream(),
      base_rng(123U),
      sampler_ptr(new sampler_t(model_, base_rng, 
                                &rstan::io::rcout, &rstan::io::rcerr)) { 
    data_stream_.close();
  }

  ~RStan() {
    delete sampler_ptr;
  }

  static void SetUpTestCase() { 
    RInside R;
  }

  static void TearDownTestCase() { }
  
  void SetUp() {
  }
  
  void TearDown() { 
  }
  

  Rcpp::List in_;
  std::fstream data_stream_;
  stan::io::dump data_context_;
  stan_model model_;
  std::fstream sample_stream;
  std::fstream diagnostic_stream;
  rng_t base_rng;
  sampler_t* sampler_ptr;
};

rstan::stan_args stan_args_factory(const std::string& str) {
  Rcpp::List list;
  std::vector<std::string> lines;
  boost::split(lines, str, boost::is_any_of("\n"));
  
  for (size_t n = 0; n < lines.size(); n++) {
    std::vector<std::string> assign;
    boost::split(assign, lines[n], boost::is_any_of("#="));
    if (assign.size() == 3) {
      std::string lhs = assign[1];
      std::string rhs = assign[2];
      boost::trim(lhs);
      boost::trim(rhs);
      
      if (lhs == "sampler_t" || lhs == "init") {
        list[lhs] = rhs;
      } else if (lhs == "seed" || lhs == "chain_id" || lhs == "iter"
                 || lhs == "warmup" || lhs == "save_warmup" || lhs == "thin" 
                 || lhs == "refresh" || lhs == "adapt_engaged" || lhs == "max_treedepth" 
                 || lhs == "append_samples") {
        int int_value;
        std::istringstream(rhs) >> int_value;
        list[lhs] = int_value;
      } else if (lhs == "stepsize" || lhs == "stepsize_jitter" || lhs == "adapt_gamma"
                 || lhs == "adapt_delta" || lhs == "adapt_kappa" 
                 || lhs == "adapt_t0") {
        double double_value;
        std::istringstream(rhs) >> double_value;
        list[lhs] = double_value;
      } else {
        ADD_FAILURE() << "don't know how to handle this line: " << lines[n];
      }
    }
  }

  return rstan::stan_args(list);
}

std::string read_stan_args(const std::string filename) {
  std::ifstream f(filename.c_str());
  if (f) {
    std::string contents;
    f.seekg(0, std::ios::end);
    contents.resize(f.tellg());
    f.seekg(0, std::ios::beg);
    f.read(&contents[0], contents.size());
    f.close();
    return(contents);
  }
  return "";
}

TEST_F(RStan, execute_sampling) {
  std::string e_stan_args_string = read_stan_args("tests/cpp/test_config/1_stan_args.txt");
  rstan::stan_args args = stan_args_factory(e_stan_args_string);

  std::stringstream ss;
  args.write_args_as_comment(ss);
  
  ASSERT_EQ(e_stan_args_string, ss.str());
}
