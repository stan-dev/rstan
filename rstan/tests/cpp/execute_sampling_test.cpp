#include <gtest/gtest.h>
#include <RInside.h>
#include <rstan/stan_fit.hpp>  
#include "blocker.cpp"
#include <sstream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <stdexcept>

typedef boost::ecuyer1988 rng_t; // (2**50 = 1T samples, 1000 chains)
typedef stan::mcmc::adapt_diag_e_nuts<stan_model, rng_t> sampler_t;

class RStan : public ::testing::Test {
public:
  RStan() 
    : in_(),
      data_stream_(std::string("tests/cpp/blocker.data.R").c_str()),
      data_context_(data_stream_),
      model_(data_context_, &std::cout),
      sample_stream(),
      diagnostic_stream(),
      base_rng1(123U),
      base_rng2(1234U),
      base_rng3(1234U),
      sampler_ptr1(new sampler_t(model_, base_rng1, 
                                 &rstan::io::rcout, &rstan::io::rcerr)),
      sampler_ptr2(new sampler_t(model_, base_rng2, 
                                 &rstan::io::rcout, &rstan::io::rcerr)),
    sampler_ptr3(new sampler_t(model_, base_rng3,
                               &rstan::io::rcout, &rstan::io::rcerr)) { 
    data_stream_.close();
  }

  ~RStan() {
    delete sampler_ptr1;
    delete sampler_ptr2;
    delete sampler_ptr3;
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
  rng_t base_rng1;
  rng_t base_rng2;
  rng_t base_rng3;
  sampler_t* sampler_ptr1;
  sampler_t* sampler_ptr2;
  sampler_t* sampler_ptr3;
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

template <class T>
T convert(const std::string& str) {
  T value;
  std::istringstream(str) >> value;
  return value;
}

template <>
std::string convert<std::string>(const std::string& str) {
  if (str[0] != '"')
    return str;
  else
    return str.substr(1, str.size()-2);
}

template <>
bool convert<bool>(const std::string& str) {
  if (str == "TRUE")
    return true;
  else if (str == "FALSE")
    return false;
  throw std::invalid_argument("not sure what \"" + str + "\" is as a boolean");
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
        list[lhs] = convert<int>(rhs);
      } else if (lhs == "stepsize" || lhs == "stepsize_jitter" || lhs == "adapt_gamma"
                 || lhs == "adapt_delta" || lhs == "adapt_kappa" 
                 || lhs == "adapt_t0") {
        list[lhs] = convert<double>(rhs);
      } else {
        ADD_FAILURE() << "don't know how to handle this line: " << lines[n];
      }
    }
  }
  return rstan::stan_args(list);
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

void parse_NumericVector(Rcpp::NumericVector& x, const std::string str) { 
  size_t start = str.find("]")+1;
  std::string values = str.substr(start);
  boost::trim(values);
  std::vector<std::string> value;
  boost::split(value, values, boost::is_any_of(" "));
  
  for (size_t m = 0; m < value.size(); m++) {
    if (value[m] != "") {
      double val;
      std::stringstream(value[m]) >> val;
      x.push_back(val);
    }
  }
}

Rcpp::List read_Rcpp_List(const std::string name, const std::vector<std::string> &lines) {
  Rcpp::List list;
  std::vector<bool> attr_lines(lines.size());
  for (size_t n = 0; n < lines.size(); n++)
    if (lines[n].substr(0, 5) == "attr(")
      attr_lines[n] = true;
    else
      attr_lines[n] = false;
  
  Rcpp::List sublist;
  std::string sublist_name;
  bool inside_sublist;
  for (size_t n = 1; n < lines.size(); n++) {
    if (n < lines.size()-1 && attr_lines[n] && attr_lines[n+1]) {
      sublist_name = lines[n].substr(lines[n].find_last_of('$') + 1);
      sublist = Rcpp::List();
      n++;
    }
    std::string lhs = lines[n].substr(lines[n].find_last_of('$') + 1);
    inside_sublist = (std::count(lines[n].begin(), lines[n].end(), '$') == 2);
    n++;
    std::string rhs = lines[n].substr(lines[n].find_last_of(']') + 1);
    boost::trim(rhs);
    
    if (sublist_name != "" && inside_sublist == false) {
      list[sublist_name] = sublist;
      sublist_name = "";
    }
    
    
    if (sublist_name == "") {
      if (lhs == "append_samples" || lhs == "test_grad")
        list[lhs] = convert<bool>(rhs);
      else if (lhs == "accept_stat__" || lhs == "stepsize__" || lhs == "treedepth__" 
               || lhs == "n_leapfrog__" || lhs == "n_divergent__") {
        Rcpp::NumericVector vec;
        while (lines[n] != "") {
          parse_NumericVector(vec, lines[n]);
          n++;
        }
        list[lhs] = vec;
      }
      else if (lhs == "chain_id" || lhs == "iter" || lhs == "refresh" 
               || lhs == "thin" || lhs == "warmup")
        list[lhs] = convert<int>(rhs);
      else if (lhs == "init" || lhs == "method" || lhs == "random_seed" || lhs == "sampler_t")
        list[lhs] = convert<std::string>(rhs);
      else if (lhs == "init_list")
        ; // this is NULL in the output. Not sure how to get that
      //list[lhs] = Rcpp::List();
      else
        list[lhs] = convert<double>(rhs);
    } else {
      if (lhs == "adapt_engaged") 
        sublist[lhs] = convert<bool>(rhs);
      else if (lhs == "adapt_init_buffer" || lhs == "adapt_t0" || lhs == "adapt_term_buffer"
               || lhs == "adapt_window" || lhs == "max_treedepth")
        sublist[lhs] = convert<int>(rhs);
      else if (lhs == "metric")
        sublist[lhs] = convert<std::string>(rhs);
      else
        sublist[lhs] = convert<double>(rhs);
    }
    
    while (n < lines.size()-1 && lines[n+1] == "")
      n++;
  }

  if (sublist_name != "")
    list[sublist_name] = sublist;
  return list;
}

std::string read_name_from_attr(const std::string str) {
  std::vector<std::string> tmp;
  boost::split(tmp, str, boost::is_any_of("\""));
  std::string name = tmp[1]; 
  boost::trim(name);
  return name;
}

Rcpp::List holder_factory(const std::string str, const rstan::stan_args args) {
  Rcpp::List holder;

  std::vector<std::string> lines;
  boost::split(lines, str, boost::is_any_of("\n"));
  
  for (int n = 0; n < lines.size(); n++) {
    if (lines[n] == "") {
      // no-op: skip line
    } else if (lines[n][0] == '$') {
      std::string name = lines[n].substr(1);
      if (name[0] == '`') {
        name = name.substr(1, name.size()-2);
      }
      n++;
      Rcpp::NumericVector x;
      while (lines[n] != "") {
        parse_NumericVector(x, lines[n]);
        n++;
      }
      holder[name] = x;
    } else if (lines[n].substr(0, 5) == "attr(") {
      std::string name = read_name_from_attr(lines[n]);
      
      if (name == "args" || name == "sampler_params") {
        std::vector<std::string> subsection;
        std::string next_name = name;
        while (next_name == name) {
          subsection.push_back(lines[n]);
          n++;
          if (n == lines.size())
            break;
          if (lines[n].substr(0, 5) == "attr(") {
            next_name = read_name_from_attr(lines[n]);
          }
        }
        n--;
        Rcpp::List x = read_Rcpp_List(name, subsection);
        holder.attr(name) = x;
      } else {
        n++;
        // NumericVector
        if (name == "inits" || name == "mean_pars") {
          Rcpp::NumericVector x;
          while (lines[n].find("]") != std::string::npos) {
            parse_NumericVector(x, lines[n]);
            n++;
          }
          n--;
          holder.attr(name) = x;
        } else if (name == "test_grad") {  // booleans
          bool x;          
          size_t start = lines[n].find("]")+1;
          std::string value = lines[n].substr(start);
          boost::trim(value);
          std::stringstream(value) >> x;
          
          holder.attr(name) = value;
        } else if (name == "mean_lp__") { // doubles
          double x;
          size_t start = lines[n].find("]")+1;
          std::string value = lines[n].substr(start);
          boost::trim(value);
          std::stringstream(value) >> x;
          
          holder.attr(name) = x;
        } else if (name == "adaptation_info") { // strings
          size_t start = lines[n].find("]")+1;
          std::string value = lines[n].substr(start);
          boost::trim(value);
          value = value.substr(1, value.size() - 2);
          boost::replace_all(value, "\\n", "\n");
          holder.attr(name) = value;
        } else {
          std::cout << "attr line: " << std::endl;
          std::cout << name << " = ";
          std::cout << "  full line: " << lines[n] << std::endl;
        }
      }
    } 
  }
  return holder;
}

void test_holder(const Rcpp::List e, const Rcpp::List x) { 
  {
    ASSERT_EQ(e.size(), x.size());
    for (size_t n = 0; n < e.size(); n++) {
      Rcpp::NumericVector e_vec = e[n];
      Rcpp::NumericVector x_vec = x[n];
      for (size_t i = 0; i < e_vec.size(); i++) {
        EXPECT_FLOAT_EQ(e_vec[i], x_vec[i]) 
          << "the " << n << "th variable, " << i << "th variable is off";
      }
    }
  }
  {
    ASSERT_TRUE(x.attr("inits") != R_NilValue);
    Rcpp::NumericVector e_inits = e.attr("inits");
    Rcpp::NumericVector x_inits = x.attr("inits");
    
    ASSERT_EQ(e_inits.size(), x_inits.size());
    for (size_t n = 0; n < e_inits.size(); n++)
      EXPECT_NEAR(e_inits[n], x_inits[n], 1e-4);
  }
  {
    ASSERT_TRUE(x.attr("mean_pars") != R_NilValue);
    Rcpp::NumericVector e_mean_pars = e.attr("mean_pars");
    Rcpp::NumericVector x_mean_pars = x.attr("mean_pars");
    
    ASSERT_EQ(e_mean_pars.size(), x_mean_pars.size());
    for (size_t n = 0; n < e_mean_pars.size(); n++)
      EXPECT_FLOAT_EQ(e_mean_pars[n], x_mean_pars[n]) 
        << "n = " << n;
  }
  {
    ASSERT_TRUE(x.attr("mean_lp__") != R_NilValue);
    double e_mean_lp__ = e.attr("mean_lp__");
    double x_mean_lp__ = x.attr("mean_lp__");
    
    EXPECT_FLOAT_EQ(e_mean_lp__, x_mean_lp__);
  }
  {
    ASSERT_TRUE(x.attr("adaptation_info") != R_NilValue);
    std::string e_adaptation_info = e.attr("adaptation_info");
    std::string x_adaptation_info = x.attr("adaptation_info");
    
    EXPECT_EQ(e_adaptation_info, x_adaptation_info);
  }
  {
    ASSERT_TRUE(x.attr("sampler_params") != R_NilValue);
    Rcpp::List e_sampler_params = e.attr("sampler_params");
    Rcpp::List x_sampler_params = x.attr("sampler_params");
    
    ASSERT_EQ(e_sampler_params.size(), x_sampler_params.size());
    
    ASSERT_TRUE(x_sampler_params["accept_stat__"] != R_NilValue);
    Rcpp::NumericVector e_accept_stat__ = e_sampler_params["accept_stat__"];
    Rcpp::NumericVector x_accept_stat__ = x_sampler_params["accept_stat__"];
    ASSERT_EQ(e_accept_stat__.size(), x_accept_stat__.size());
    for (size_t n = 0; n < e_accept_stat__.size(); n++)
      EXPECT_NEAR(e_accept_stat__[n], x_accept_stat__[n], 1e-7);

    ASSERT_TRUE(x_sampler_params["stepsize__"] != R_NilValue);
    Rcpp::NumericVector e_stepsize__ = e_sampler_params["stepsize__"];
    Rcpp::NumericVector x_stepsize__ = x_sampler_params["stepsize__"];
    ASSERT_EQ(e_stepsize__.size(), x_stepsize__.size());
    for (size_t n = 0; n < e_stepsize__.size(); n++)
      EXPECT_FLOAT_EQ(e_stepsize__[n], x_stepsize__[n]);

    ASSERT_TRUE(x_sampler_params["treedepth__"] != R_NilValue);
    Rcpp::NumericVector e_treedepth__ = e_sampler_params["treedepth__"];
    Rcpp::NumericVector x_treedepth__ = x_sampler_params["treedepth__"];
    ASSERT_EQ(e_treedepth__.size(), x_treedepth__.size());
    for (size_t n = 0; n < e_treedepth__.size(); n++)
      EXPECT_FLOAT_EQ(e_treedepth__[n], x_treedepth__[n]);

    ASSERT_TRUE(x_sampler_params["n_leapfrog__"] != R_NilValue);
    Rcpp::NumericVector e_n_leapfrog__ = e_sampler_params["n_leapfrog__"];
    Rcpp::NumericVector x_n_leapfrog__ = x_sampler_params["n_leapfrog__"];
    ASSERT_EQ(e_n_leapfrog__.size(), x_n_leapfrog__.size());
    for (size_t n = 0; n < e_n_leapfrog__.size(); n++)
      EXPECT_FLOAT_EQ(e_n_leapfrog__[n], x_n_leapfrog__[n]);

    ASSERT_TRUE(x_sampler_params["n_divergent__"] != R_NilValue);
    Rcpp::NumericVector e_n_divergent__ = e_sampler_params["n_divergent__"];
    Rcpp::NumericVector x_n_divergent__ = x_sampler_params["n_divergent__"];
    ASSERT_EQ(e_n_divergent__.size(), x_n_divergent__.size());
    for (size_t n = 0; n < e_n_divergent__.size(); n++)
      EXPECT_FLOAT_EQ(e_n_divergent__[n], x_n_divergent__[n]);
  }
}



TEST_F(RStan, execute_sampling_1) {

  std::stringstream ss;

  std::string e_args_string = read_file("tests/cpp/test_config/1_input_stan_args.txt");
  rstan::stan_args args = stan_args_factory(e_args_string);
  args.write_args_as_comment(ss);
  ASSERT_EQ(e_args_string, ss.str());
  ss.str("");
  
  stan::mcmc::sample s(Eigen::VectorXd::Zero(model_.num_params_r()), 0, 0);
  
  std::vector<size_t> qoi_idx 
    = vector_factory<size_t>(read_file("tests/cpp/test_config/1_input_qoi_idx.txt"));
  std::vector<double> initv 
    = vector_factory<double>(read_file("tests/cpp/test_config/1_input_initv.txt"));
  std::vector<std::string> fnames_oi 
    = vector_factory<std::string>(read_file("tests/cpp/test_config/1_input_fnames_oi.txt"));
  
  std::string e_holder_string = read_file("tests/cpp/test_config/1_output_holder.txt");
  Rcpp::List e_holder 
    = holder_factory(read_file("tests/cpp/test_config/1_output_holder.txt"),
                     args);

  Rcpp::List holder;
  
  rstan::init_nuts<sampler_t>(sampler_ptr1, args);
  Eigen::VectorXd tmp = Eigen::VectorXd::Zero(model_.num_params_r());
  rstan::init_windowed_adapt<sampler_t>(sampler_ptr1, args, s.cont_params());
  
  
  std::fstream sample_stream, diagnostic_stream;
  rstan::execute_sampling(args, model_, holder,
                          sampler_ptr1,
                          s,
                          qoi_idx,
                          initv,
                          sample_stream,
                          diagnostic_stream,
                          fnames_oi,
                          base_rng1);
  
  test_holder(e_holder, holder);
}


TEST_F(RStan, execute_sampling_2) {
  std::stringstream ss;

  std::string e_args_string = read_file("tests/cpp/test_config/2_input_stan_args.txt");
  rstan::stan_args args = stan_args_factory(e_args_string);
  args.write_args_as_comment(ss);
  ASSERT_EQ(e_args_string, ss.str());
  ss.str("");
  
  stan::mcmc::sample s(Eigen::VectorXd::Zero(model_.num_params_r()), 0, 0);
  
  std::vector<size_t> qoi_idx 
    = vector_factory<size_t>(read_file("tests/cpp/test_config/2_input_qoi_idx.txt"));
  std::vector<double> initv 
    = vector_factory<double>(read_file("tests/cpp/test_config/2_input_initv.txt"));
  std::vector<std::string> fnames_oi 
    = vector_factory<std::string>(read_file("tests/cpp/test_config/2_input_fnames_oi.txt"));
  
  std::string e_holder_string = read_file("tests/cpp/test_config/2_output_holder.txt");
  Rcpp::List e_holder 
    = holder_factory(read_file("tests/cpp/test_config/2_output_holder.txt"),
                     args);

  Rcpp::List holder;
  
  rstan::init_nuts<sampler_t>(sampler_ptr2, args);
  Eigen::VectorXd tmp = Eigen::VectorXd::Zero(model_.num_params_r());
  rstan::init_windowed_adapt<sampler_t>(sampler_ptr2, args, s.cont_params());
  
  
  std::fstream sample_stream, diagnostic_stream;
  rstan::execute_sampling(args, model_, holder,
                          sampler_ptr2,
                          s,
                          qoi_idx,
                          initv,
                          sample_stream,
                          diagnostic_stream,
                          fnames_oi,
                          base_rng2);
  
  
  
  test_holder(e_holder, holder);
}

TEST_F(RStan, execute_sampling_3) {
  std::stringstream ss;

  std::string e_args_string = read_file("tests/cpp/test_config/3_input_stan_args.txt");
  rstan::stan_args args = stan_args_factory(e_args_string);
  args.write_args_as_comment(ss);
  ASSERT_EQ(e_args_string, ss.str());
  ss.str("");

  double r = args.get_init_radius();
  boost::random::uniform_real_distribution<double> 
    init_range_distribution(-r, r);
  boost::variate_generator<rng_t&, boost::random::uniform_real_distribution<double> >
    init_rng(base_rng3,init_range_distribution);
  Eigen::VectorXd cont_params(model_.num_params_r());
  for (size_t i = 0; i < cont_params.size(); ++i)
    cont_params[i] = init_rng();
  
  stan::mcmc::sample s(cont_params, 0, 0);
  
  
  std::vector<size_t> qoi_idx 
    = vector_factory<size_t>(read_file("tests/cpp/test_config/3_input_qoi_idx.txt"));
  std::vector<double> initv 
    = vector_factory<double>(read_file("tests/cpp/test_config/3_input_initv.txt"));
  std::vector<std::string> fnames_oi 
    = vector_factory<std::string>(read_file("tests/cpp/test_config/3_input_fnames_oi.txt"));
  
  std::string e_holder_string = read_file("tests/cpp/test_config/3_output_holder.txt");
  Rcpp::List e_holder 
    = holder_factory(read_file("tests/cpp/test_config/3_output_holder.txt"),
                     args);

  Rcpp::List holder;
  
  rstan::init_nuts<sampler_t>(sampler_ptr2, args);
  Eigen::VectorXd tmp = Eigen::VectorXd::Zero(model_.num_params_r());
  rstan::init_windowed_adapt<sampler_t>(sampler_ptr2, args, s.cont_params());
  
  
  std::fstream sample_stream, diagnostic_stream;
  rstan::execute_sampling(args, model_, holder,
                          sampler_ptr2,
                          s,
                          qoi_idx,
                          initv,
                          sample_stream,
                          diagnostic_stream,
                          fnames_oi,
                          base_rng2);
  
  
  
  test_holder(e_holder, holder);
}
