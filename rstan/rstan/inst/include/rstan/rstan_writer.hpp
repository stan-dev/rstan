#ifndef RSTAN__RSTAN_WRITER_HPP
#define RSTAN__RSTAN_WRITER_HPP

#include <Rcpp.h>
#include <stan/callbacks/writer.hpp>
#include <stan/callbacks/stream_writer.hpp>
#include <rstan/filtered_values.hpp>
#include <rstan/sum_values.hpp>

namespace rstan {

  class rstan_sample_writer : public stan::callbacks::writer {
  public:
    stan::callbacks::stream_writer csv_;
    filtered_values<Rcpp::NumericVector> values_;
    filtered_values<Rcpp::NumericVector> sampler_values_;
    sum_values sum_;

    rstan_sample_writer(stan::callbacks::stream_writer csv,
                        filtered_values<Rcpp::NumericVector> values,
                        filtered_values<Rcpp::NumericVector> sampler_values,
                        sum_values sum)
      : csv_(csv), values_(values), sampler_values_(sampler_values), sum_(sum) { }

    /**
     * Writes a key, value pair.
     *
     * @param[in] key A string
     * @param[in] value A double value
     */
    void operator()(const std::string& key, double value) {
      csv_(key, value);
      values_(key, value);
      sampler_values_(key, value);
      sum_(key, value);
    }

    /**
     * Writes a key, value pair.
     *
     * @param[in] key A string
     * @param[in] value An integer value
     */
    void operator()(const std::string& key, int value) {
      csv_(key, value);
      values_(key, value);
      sampler_values_(key, value);
      sum_(key, value);
    }

    /**
     * Writes a key, value pair.
     *
     * @param[in] key A string
     * @param[in] value A string
     */
    void operator()(const std::string& key, const std::string& value) {
      csv_(key, value);
      values_(key, value);
      sampler_values_(key, value);
      sum_(key, value);
    }

    /**
     * Writes a key, value pair.
     *
     * @param[in] key A string
     * @param[in] values A double array, typically used with
     *   contiguous Eigen vectors
     * @param[in] n_values Length of the array
     */
    void operator()(const std::string& key, const double* values,
                    int n_values)  {
      csv_(key, values, n_values);
      values_(key, values, n_values);
      sampler_values_(key, values, n_values);
      sum_(key, values, n_values);
    }

    /**
     * Writes a key, value pair.
     *
     * @param[in] key A string
     * @param[in] values A double array assumed to represent a 2d
     *   matrix stored in column major order, typically used with
     *   contiguous Eigen matrices
     * @param[in] n_rows Rows
     * @param[in] n_cols Columns
     */
    void operator()(const std::string& key, const double* values,
                    int n_rows, int n_cols) {
      csv_(key, values, n_rows, n_cols);
      values_(key, values, n_rows, n_cols);
      sampler_values_(key, values, n_rows, n_cols);
      sum_(key, values, n_rows, n_cols);
    }

    /**
     * Writes a set of names.
     *
     * @param[in] names Names in a std::vector
     */
    void operator()(const std::vector<std::string>& names) {
      csv_(names);
      values_(names);
      sampler_values_(names);
      sum_(names);
    }

    /**
     * Writes a set of values.
     *
     * @param[in] state Values in a std::vector
     */
    void operator()(const std::vector<double>& state) {
      csv_(state);
      values_(state);
      sampler_values_(state);
      sum_(state);
    }

    /**
     * Writes a string.
     *
     * @param[in] message A string
     */
    void operator()(const std::string& message) {
      csv_(message);
      values_(message);
      sampler_values_(message);
      sum_(message);
    }

    /**
     * Writes blank input.
     */
    void operator()() {
      csv_();
      values_();
      sampler_values_();
      sum_();
    }
  };

  /**
     @param      N
     @param      M  number of iterations to be saved
     @param      warmup number of warmup iterations to be saved
  */
  rstan_sample_writer*
  sample_writer_factory(std::ostream *o, const std::string& prefix,
                        size_t N_sample_names, size_t N_sampler_names,
                        size_t N_constrained_param_names,
                        size_t N_iter_save, size_t warmup,
                        size_t offset,
                        const std::vector<size_t>& qoi_idx) {
    size_t N = N_sample_names + N_sampler_names + N_constrained_param_names;
    
    std::vector<size_t> filter(qoi_idx);
    std::vector<size_t> lp;
    for (size_t n = 0; n < filter.size(); n++)
      if (filter[n] >= N)
        lp.push_back(n);
    for (size_t n = 0; n < filter.size(); n++)
      filter[n] += offset;
    for (size_t n = 0; n < lp.size(); n++)
      filter[lp[n]] = 0;

    std::vector<size_t> filter_sampler_values(offset);
    for (size_t n = 0; n < offset; n++)
      filter_sampler_values[n] = n;

    stan::callbacks::stream_writer csv(*o, prefix);
    filtered_values<Rcpp::NumericVector> values(N, N_iter_save, filter);
    filtered_values<Rcpp::NumericVector> sampler_values(N, N_iter_save, filter_sampler_values);
    sum_values sum(N, warmup);

    return new rstan_sample_writer(csv, values, sampler_values, sum);
  }

}
#endif
