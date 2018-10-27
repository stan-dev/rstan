#ifndef RSTAN__RSTAN_WRITER_HPP
#define RSTAN__RSTAN_WRITER_HPP

#include <Rcpp.h>
#include <stan/callbacks/writer.hpp>
#include <stan/callbacks/stream_writer.hpp>
#include <rstan/comment_writer.hpp>
#include <rstan/filtered_values.hpp>
#include <rstan/sum_values.hpp>

namespace rstan {

  class rstan_sample_writer : public stan::callbacks::writer {
  public:
    stan::callbacks::stream_writer csv_;
    comment_writer comment_writer_;
    filtered_values<Rcpp::NumericVector> values_;
    filtered_values<Rcpp::NumericVector> sampler_values_;
    sum_values sum_;

    rstan_sample_writer(stan::callbacks::stream_writer csv,
                        comment_writer comment_writer,
                        filtered_values<Rcpp::NumericVector> values,
                        filtered_values<Rcpp::NumericVector> sampler_values,
                        sum_values sum)
      : csv_(csv), comment_writer_(comment_writer),
        values_(values), sampler_values_(sampler_values), sum_(sum) { }

    /**
     * Writes a set of names.
     *
     * @param[in] names Names in a std::vector
     */
    void operator()(const std::vector<std::string>& names) {
      csv_(names);
      comment_writer_(names);
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
      comment_writer_(state);
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
      comment_writer_(message);
      values_(message);
      sampler_values_(message);
      sum_(message);
    }

    /**
     * Writes blank input.
     */
    void operator()() {
      csv_();
      comment_writer_();
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
  inline
  rstan_sample_writer*
  sample_writer_factory(std::ostream *csv_fstream,
                        std::ostream& comment_stream,
                        const std::string& prefix,
                        size_t N_sample_names, size_t N_sampler_names,
                        size_t N_constrained_param_names,
                        size_t N_iter_save, size_t warmup,
                        const std::vector<size_t>& qoi_idx) {
    size_t N = N_sample_names + N_sampler_names + N_constrained_param_names;
    size_t offset = N_sample_names + N_sampler_names;

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

    stan::callbacks::stream_writer csv(*csv_fstream, prefix);
    comment_writer comments(comment_stream, prefix);
    filtered_values<Rcpp::NumericVector> values(N, N_iter_save, filter);
    filtered_values<Rcpp::NumericVector> sampler_values(N, N_iter_save, filter_sampler_values);
    sum_values sum(N, warmup);

    return new rstan_sample_writer(csv, comments, values, sampler_values, sum);
  }

}
#endif
