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

std::string read_file(const std::string filename) {
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

template <class T>
T convert(const std::string& str) {
  T value;
  std::istringstream(str) >> value;
  return value;
}

template <>
std::string convert(const std::string& str) {
  return str;
}

template <class T>
std::vector<T> vector_factory(const std::string str) {
  std::vector<T> x;
  std::vector<std::string> lines;
  boost::split(lines, str, boost::is_any_of("\n"));
  if (lines.size() == 3) {
    size_t size;
    std::vector<std::string> line;
    boost::split(line, lines[0], boost::is_any_of("()"));
    if (line.size() < 2) {
      ADD_FAILURE() << "can't find the size: " << str;
      return x;
    }
    std::istringstream(line[1]) >> size;
    x.resize(size);
    

    line.clear();
    boost::split(line, lines[1], boost::is_any_of("{},"));
    if (line.size() != size + 2) {
      ADD_FAILURE() << "can't read the correct number of elements: " << str;
      return x;
    }
    
    for (size_t n = 0; n < size; n++) {
      x[n] = convert<T>(line[n+1]);
    }
  } else {
    ADD_FAILURE() << "number of lines (" << lines.size() << ") is not 3: "
                  << str << std::endl;
  }
  return x;
}


TEST_F(RStan, execute_sampling_1) {
  std::stringstream ss;

  std::string e_args_string = read_file("tests/cpp/test_config/1_input_stan_args.txt");
  rstan::stan_args args = stan_args_factory(e_args_string);
  args.write_args_as_comment(ss);
  ASSERT_EQ(e_args_string, ss.str());
  ss.str("");
  
  stan::mcmc::sample s(Eigen::VectorXd(model_.num_params_r()), 0, 0);
  
  std::vector<size_t> qoi_idx 
    = vector_factory<size_t>(read_file("tests/cpp/test_config/1_input_qoi_idx.txt"));
  std::vector<double> initv 
    = vector_factory<double>(read_file("tests/cpp/test_config/1_input_initv.txt"));
  std::vector<std::string> fnames_oi 
    = vector_factory<std::string>(read_file("tests/cpp/test_config/1_input_fnames_oi.txt"));
  
  
}
