# RStan — the R interface to Stan

*Stan Development Team*

RStan is the R interface to the [Stan](https://mc-stan.org/) C++
package. The RStan interface (rstan R package) provides:

- Full Bayesian inference using the No-U-Turn sampler (NUTS), a variant
  of Hamiltonian Monte Carlo (HMC)

- Approximate Bayesian inference using automatic differentiation
  variational inference (ADVI)

- Penalized maximum likelihood estimation using L-BFGS optimization

For documentation on Stan itself, including the manual and user guide
for the modeling language, case studies and worked examples, and other
tutorial information visit the Users section of the Stan website:

- [mc-stan.org/users/documentation](https://mc-stan.org/users/documentation/)

## Other R packages from the Stan Development Team

Various related R packages are also available from the Stan Development
Team including these and more:

|  |  |  |  |
|----|----|----|----|
| **Package** | **Description** | **Doc** | **Website** |
| bayesplot | ggplot-based plotting of parameter estimates, diagnostics, and posterior predictive checks. | [bayesplot-package](https://mc-stan.org/bayesplot/reference/bayesplot-package.html) | [mc-stan.org/bayesplot](https://mc-stan.org/bayesplot/) |
| shinystan | Interactive GUI for exploring MCMC output. | [shinystan-package](https://mc-stan.org/shinystan/reference/shinystan-package.html) | [mc-stan.org/shinystan](https://mc-stan.org/shinystan/) |
| loo | Out-of-sample predictive performance estimates and model comparison. | [loo-package](https://mc-stan.org/loo/reference/loo-package.html) | [mc-stan.org/loo](https://mc-stan.org/loo/) |
| rstanarm | R formula interface for applied regression modeling. | rstanarm-package | [mc-stan.org/rstanarm](https://mc-stan.org/rstanarm/) |
| rstantools | Tools for developers of R packages interfacing with Stan. | [rstantools-package](https://mc-stan.org/rstantools/reference/rstantools-package.html) | [mc-stan.org/rstantools](https://mc-stan.org/rstantools/) |

## Author

|                                   |                                    |
|-----------------------------------|------------------------------------|
| Jonah Gabry (author)              | \<jonah.sol.gabry@columbia.edu\>   |
| Ben Goodrich (maintainer, author) | \<benjamin.goodrich@columbia.edu\> |
| Jiqiang Guo (author)              | \<guojq28@gmail.com\>              |

There are also many other important contributors to RStan
([github.com/rstan](https://github.com/stan-dev/rstan)). Please use
'Stan Development Team' whenever citing the R interface to Stan. A
BibTex entry is available from <https://mc-stan.org/rstan/authors> or
`citation("rstan")`.

## See also

- The RStan vignettes: <https://mc-stan.org/rstan/articles/>.

- [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md) for details
  on fitting models and
  [`stanfit`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md)
  for information on the fitted model objects.

- The [`lookup`](https://mc-stan.org/rstan/dev/reference/lookup.md) for
  finding a function in the Stan language that corresponds to a R
  function or name.

- <https://github.com/stan-dev/rstan/issues/> to submit a bug report or
  feature request.

- <https://discourse.mc-stan.org> to ask a question on the Stan Forums.

## Examples

``` r
# \dontrun{

stanmodelcode <- "
data {
  int<lower=0> N;
  array[N] real y;
}

parameters {
  real mu;
}

model {
  target += normal_lpdf(mu | 0, 10);
  target += normal_lpdf(y  | mu, 1);
}
"

y <- rnorm(20)
dat <- list(N = 20, y = y);
fit <- stan(model_code = stanmodelcode, model_name = "example",
            data = dat, iter = 2012, chains = 3, verbose = TRUE,
            sample_file = file.path(tempdir(), 'norm.csv'))
#> 
#> TRANSLATING MODEL 'example' FROM Stan CODE TO C++ CODE NOW.
#> OS: x86_64, linux-gnu; rstan: 2.39.0.9000; Rcpp: 1.1.2; inline: 0.3.21 
#>  >> setting environment variables: 
#> PKG_LIBS =  '/home/runner/work/_temp/Library/rstan/lib//libStanServices.a' -L'/home/runner/work/_temp/Library/StanHeaders/lib/' -lStanHeaders -L'/home/runner/work/_temp/Library/RcppParallel/lib/' -ltbb 
#> PKG_CPPFLAGS =   -I"/home/runner/work/_temp/Library/Rcpp/include/"  -I"/home/runner/work/_temp/Library/RcppEigen/include/"  -I"/home/runner/work/_temp/Library/RcppEigen/include/unsupported"  -I"/home/runner/work/_temp/Library/BH/include" -I"/home/runner/work/_temp/Library/StanHeaders/include/src/"  -I"/home/runner/work/_temp/Library/StanHeaders/include/"  -I"/home/runner/work/_temp/Library/RcppParallel/include/" -DRCPP_PARALLEL_USE_TBB=1 -DTBB_INTERFACE_NEW -I/home/runner/work/_temp/Library/RcppParallel/include -I"/home/runner/work/_temp/Library/rstan/include" -DEIGEN_NO_DEBUG  -DBOOST_DISABLE_ASSERTS  -DBOOST_PENDING_INTEGER_LOG2_HPP  -DSTAN_THREADS  -DUSE_STANC3 -DSTRICT_R_HEADERS  -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION  -D_HAS_AUTO_PTR_ETC=0  -include '/home/runner/work/_temp/Library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp'  -D_REENTRANT -DRCPP_PARALLEL_USE_TBB=1
#>  >> Program source :
#> 
#>    1 : 
#>    2 : // includes from the plugin
#>    3 : // [[Rcpp::plugins(cpp17)]]
#>    4 : 
#>    5 : 
#>    6 : // user includes
#>    7 : #include <Rcpp.h>
#>    8 : using namespace Rcpp;
#>    9 : #ifndef MODELS_HPP
#>   10 : #define MODELS_HPP
#>   11 : #define STAN__SERVICES__COMMAND_HPP
#>   12 : #include <rstan/rstaninc.hpp>
#>   13 : #ifndef USE_STANC3
#>   14 : #define USE_STANC3
#>   15 : #endif
#>   16 : // Code generated by stanc 5496e41
#>   17 : #include <stan/model/model_header.hpp>
#>   18 : namespace model1c8428b78e6e_example_namespace {
#>   19 : using stan::model::model_base_crtp;
#>   20 : using namespace stan::math;
#>   21 : stan::math::profile_map profiles__;
#>   22 : static constexpr std::array<const char*, 7> locations_array__ =
#>   23 :   {" (found before start of program)",
#>   24 :   " (in 'example', line 6, column 2 to column 10)",
#>   25 :   " (in 'example', line 9, column 2 to column 36)",
#>   26 :   " (in 'example', line 10, column 2 to column 36)",
#>   27 :   " (in 'example', line 2, column 2 to column 17)",
#>   28 :   " (in 'example', line 3, column 8 to column 9)",
#>   29 :   " (in 'example', line 3, column 2 to column 18)"};
#>   30 : class model1c8428b78e6e_example final : public model_base_crtp<model1c8428b78e6e_example> {
#>   31 : private:
#>   32 :   int N;
#>   33 :   std::vector<double> y;
#>   34 : public:
#>   35 :   ~model1c8428b78e6e_example() {}
#>   36 :   model1c8428b78e6e_example(stan::io::var_context& context__, unsigned int
#>   37 :                             random_seed__ = 0, std::ostream*
#>   38 :                             pstream__ = nullptr) : model_base_crtp(0) {
#>   39 :     int current_statement__ = 0;
#>   40 :     // suppress unused var warning
#>   41 :     (void) current_statement__;
#>   42 :     using local_scalar_t__ = double;
#>   43 :     auto base_rng__ = stan::services::util::create_rng(random_seed__, 0);
#>   44 :     // suppress unused var warning
#>   45 :     (void) base_rng__;
#>   46 :     static constexpr const char* function__ =
#>   47 :       "model1c8428b78e6e_example_namespace::model1c8428b78e6e_example";
#>   48 :     // suppress unused var warning
#>   49 :     (void) function__;
#>   50 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>   51 :     // suppress unused var warning
#>   52 :     (void) DUMMY_VAR__;
#>   53 :     try {
#>   54 :       current_statement__ = 4;
#>   55 :       context__.validate_dims("data initialization", "N", "int",
#>   56 :         std::vector<size_t>{});
#>   57 :       N = std::numeric_limits<int>::min();
#>   58 :       current_statement__ = 4;
#>   59 :       N = context__.vals_i("N")[(1 - 1)];
#>   60 :       current_statement__ = 4;
#>   61 :       stan::math::check_greater_or_equal(function__, "N", N, 0);
#>   62 :       current_statement__ = 5;
#>   63 :       stan::math::validate_non_negative_index("y", "N", N);
#>   64 :       current_statement__ = 6;
#>   65 :       context__.validate_dims("data initialization", "y", "double",
#>   66 :         std::vector<size_t>{static_cast<size_t>(N)});
#>   67 :       y = std::vector<double>(N, std::numeric_limits<double>::quiet_NaN());
#>   68 :       current_statement__ = 6;
#>   69 :       y = context__.vals_r("y");
#>   70 :     } catch (const std::exception& e) {
#>   71 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>   72 :     }
#>   73 :     num_params_r__ = 1;
#>   74 :   }
#>   75 :   inline std::string model_name() const final {
#>   76 :     return "model1c8428b78e6e_example";
#>   77 :   }
#>   78 :   inline std::vector<std::string> model_compile_info() const noexcept {
#>   79 :     return std::vector<std::string>{"stanc_version = stanc3 5496e41",
#>   80 :              "stancflags = --"};
#>   81 :   }
#>   82 :   // Base log prob
#>   83 :   template <bool propto__, bool jacobian__, typename VecR, typename VecI,
#>   84 :             stan::require_vector_like_t<VecR>* = nullptr,
#>   85 :             stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr,
#>   86 :             stan::require_not_st_var<VecR>* = nullptr>
#>   87 :   inline stan::scalar_type_t<VecR>
#>   88 :   log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
#>   89 :                 pstream__ = nullptr) const {
#>   90 :     using T__ = stan::scalar_type_t<VecR>;
#>   91 :     using local_scalar_t__ = T__;
#>   92 :     T__ lp__(0.0);
#>   93 :     stan::math::accumulator<T__> lp_accum__;
#>   94 :     stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
#>   95 :     int current_statement__ = 0;
#>   96 :     // suppress unused var warning
#>   97 :     (void) current_statement__;
#>   98 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>   99 :     // suppress unused var warning
#>  100 :     (void) DUMMY_VAR__;
#>  101 :     static constexpr const char* function__ =
#>  102 :       "model1c8428b78e6e_example_namespace::log_prob";
#>  103 :     // suppress unused var warning
#>  104 :     (void) function__;
#>  105 :     try {
#>  106 :       current_statement__ = 1;
#>  107 :       auto mu = in__.template read<local_scalar_t__>();
#>  108 :       {
#>  109 :         current_statement__ = 2;
#>  110 :         lp_accum__.add(stan::math::normal_lpdf<false>(mu,
#>  111 :                          static_cast<double>(0), static_cast<double>(10)));
#>  112 :         current_statement__ = 3;
#>  113 :         lp_accum__.add(stan::math::normal_lpdf<false>(y, mu,
#>  114 :                          static_cast<double>(1)));
#>  115 :       }
#>  116 :     } catch (const std::exception& e) {
#>  117 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  118 :     }
#>  119 :     lp_accum__.add(lp__);
#>  120 :     return lp_accum__.sum();
#>  121 :   }
#>  122 :   // Reverse mode autodiff log prob
#>  123 :   template <bool propto__, bool jacobian__, typename VecR, typename VecI,
#>  124 :             stan::require_vector_like_t<VecR>* = nullptr,
#>  125 :             stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr,
#>  126 :             stan::require_st_var<VecR>* = nullptr>
#>  127 :   inline stan::scalar_type_t<VecR>
#>  128 :   log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
#>  129 :                 pstream__ = nullptr) const {
#>  130 :     using T__ = stan::scalar_type_t<VecR>;
#>  131 :     using local_scalar_t__ = T__;
#>  132 :     T__ lp__(0.0);
#>  133 :     stan::math::accumulator<T__> lp_accum__;
#>  134 :     stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
#>  135 :     int current_statement__ = 0;
#>  136 :     // suppress unused var warning
#>  137 :     (void) current_statement__;
#>  138 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>  139 :     // suppress unused var warning
#>  140 :     (void) DUMMY_VAR__;
#>  141 :     static constexpr const char* function__ =
#>  142 :       "model1c8428b78e6e_example_namespace::log_prob";
#>  143 :     // suppress unused var warning
#>  144 :     (void) function__;
#>  145 :     try {
#>  146 :       current_statement__ = 1;
#>  147 :       auto mu = in__.template read<local_scalar_t__>();
#>  148 :       {
#>  149 :         current_statement__ = 2;
#>  150 :         lp_accum__.add(stan::math::normal_lpdf<false>(mu,
#>  151 :                          static_cast<double>(0), static_cast<double>(10)));
#>  152 :         current_statement__ = 3;
#>  153 :         lp_accum__.add(stan::math::normal_lpdf<false>(y, mu,
#>  154 :                          static_cast<double>(1)));
#>  155 :       }
#>  156 :     } catch (const std::exception& e) {
#>  157 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  158 :     }
#>  159 :     lp_accum__.add(lp__);
#>  160 :     return lp_accum__.sum();
#>  161 :   }
#>  162 :   template <typename RNG, typename VecR, typename VecI, typename VecVar,
#>  163 :             stan::require_vector_like_vt<std::is_floating_point,
#>  164 :             VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
#>  165 :             VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
#>  166 :             VecVar>* = nullptr>
#>  167 :   inline void
#>  168 :   write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
#>  169 :                    VecVar& vars__, const bool
#>  170 :                    emit_transformed_parameters__ = true, const bool
#>  171 :                    emit_generated_quantities__ = true, std::ostream*
#>  172 :                    pstream__ = nullptr) const {
#>  173 :     using local_scalar_t__ = double;
#>  174 :     stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
#>  175 :     stan::io::serializer<local_scalar_t__> out__(vars__);
#>  176 :     static constexpr bool propto__ = true;
#>  177 :     // suppress unused var warning
#>  178 :     (void) propto__;
#>  179 :     double lp__ = 0.0;
#>  180 :     // suppress unused var warning
#>  181 :     (void) lp__;
#>  182 :     int current_statement__ = 0;
#>  183 :     // suppress unused var warning
#>  184 :     (void) current_statement__;
#>  185 :     stan::math::accumulator<double> lp_accum__;
#>  186 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>  187 :     // suppress unused var warning
#>  188 :     (void) DUMMY_VAR__;
#>  189 :     constexpr bool jacobian__ = false;
#>  190 :     // suppress unused var warning
#>  191 :     (void) jacobian__;
#>  192 :     static constexpr const char* function__ =
#>  193 :       "model1c8428b78e6e_example_namespace::write_array";
#>  194 :     // suppress unused var warning
#>  195 :     (void) function__;
#>  196 :     try {
#>  197 :       current_statement__ = 1;
#>  198 :       auto mu = in__.template read<local_scalar_t__>();
#>  199 :       out__.write(mu);
#>  200 :       if (stan::math::logical_negation(
#>  201 :             (stan::math::primitive_value(emit_transformed_parameters__) ||
#>  202 :             stan::math::primitive_value(emit_generated_quantities__)))) {
#>  203 :         return ;
#>  204 :       }
#>  205 :       if (stan::math::logical_negation(emit_generated_quantities__)) {
#>  206 :         return ;
#>  207 :       }
#>  208 :     } catch (const std::exception& e) {
#>  209 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  210 :     }
#>  211 :   }
#>  212 :   template <typename VecVar, typename VecI,
#>  213 :             stan::require_vector_t<VecVar>* = nullptr,
#>  214 :             stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
#>  215 :   inline void
#>  216 :   unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
#>  217 :                          VecVar& vars__, std::ostream* pstream__ = nullptr) const {
#>  218 :     using local_scalar_t__ = double;
#>  219 :     stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
#>  220 :     stan::io::serializer<local_scalar_t__> out__(vars__);
#>  221 :     int current_statement__ = 0;
#>  222 :     // suppress unused var warning
#>  223 :     (void) current_statement__;
#>  224 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>  225 :     // suppress unused var warning
#>  226 :     (void) DUMMY_VAR__;
#>  227 :     try {
#>  228 :       local_scalar_t__ mu = DUMMY_VAR__;
#>  229 :       current_statement__ = 1;
#>  230 :       mu = in__.read<local_scalar_t__>();
#>  231 :       out__.write(mu);
#>  232 :     } catch (const std::exception& e) {
#>  233 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  234 :     }
#>  235 :   }
#>  236 :   template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
#>  237 :   inline void
#>  238 :   transform_inits_impl(const stan::io::var_context& context__, VecVar&
#>  239 :                        vars__, std::ostream* pstream__ = nullptr) const {
#>  240 :     using local_scalar_t__ = double;
#>  241 :     stan::io::serializer<local_scalar_t__> out__(vars__);
#>  242 :     int current_statement__ = 0;
#>  243 :     // suppress unused var warning
#>  244 :     (void) current_statement__;
#>  245 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>  246 :     // suppress unused var warning
#>  247 :     (void) DUMMY_VAR__;
#>  248 :     try {
#>  249 :       current_statement__ = 1;
#>  250 :       context__.validate_dims("parameter initialization", "mu", "double",
#>  251 :         std::vector<size_t>{});
#>  252 :       local_scalar_t__ mu = DUMMY_VAR__;
#>  253 :       current_statement__ = 1;
#>  254 :       mu = context__.vals_r("mu")[(1 - 1)];
#>  255 :       out__.write(mu);
#>  256 :     } catch (const std::exception& e) {
#>  257 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  258 :     }
#>  259 :   }
#>  260 :   inline void
#>  261 :   get_param_names(std::vector<std::string>& names__, const bool
#>  262 :                   emit_transformed_parameters__ = true, const bool
#>  263 :                   emit_generated_quantities__ = true) const {
#>  264 :     names__ = std::vector<std::string>{"mu"};
#>  265 :     if (emit_transformed_parameters__) {}
#>  266 :     if (emit_generated_quantities__) {}
#>  267 :   }
#>  268 :   inline void
#>  269 :   get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
#>  270 :            emit_transformed_parameters__ = true, const bool
#>  271 :            emit_generated_quantities__ = true) const {
#>  272 :     dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{}};
#>  273 :     if (emit_transformed_parameters__) {}
#>  274 :     if (emit_generated_quantities__) {}
#>  275 :   }
#>  276 :   inline void
#>  277 :   constrained_param_names(std::vector<std::string>& param_names__, bool
#>  278 :                           emit_transformed_parameters__ = true, bool
#>  279 :                           emit_generated_quantities__ = true) const final {
#>  280 :     param_names__.emplace_back(std::string() + "mu");
#>  281 :     if (emit_transformed_parameters__) {}
#>  282 :     if (emit_generated_quantities__) {}
#>  283 :   }
#>  284 :   inline void
#>  285 :   unconstrained_param_names(std::vector<std::string>& param_names__, bool
#>  286 :                             emit_transformed_parameters__ = true, bool
#>  287 :                             emit_generated_quantities__ = true) const final {
#>  288 :     param_names__.emplace_back(std::string() + "mu");
#>  289 :     if (emit_transformed_parameters__) {}
#>  290 :     if (emit_generated_quantities__) {}
#>  291 :   }
#>  292 :   inline std::string get_constrained_sizedtypes() const {
#>  293 :     return std::string("[{\"name\":\"mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]");
#>  294 :   }
#>  295 :   inline std::string get_unconstrained_sizedtypes() const {
#>  296 :     return std::string("[{\"name\":\"mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]");
#>  297 :   }
#>  298 :   // Begin method overload boilerplate
#>  299 :   template <typename RNG> inline void
#>  300 :   write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
#>  301 :               Eigen::Matrix<double,-1,1>& vars, const bool
#>  302 :               emit_transformed_parameters = true, const bool
#>  303 :               emit_generated_quantities = true, std::ostream*
#>  304 :               pstream = nullptr) const {
#>  305 :     const size_t num_params__ = 1;
#>  306 :     const size_t num_transformed = emit_transformed_parameters * (0U);
#>  307 :     const size_t num_gen_quantities = emit_generated_quantities * (0U);
#>  308 :     const size_t num_to_write = num_params__ + num_transformed +
#>  309 :       num_gen_quantities;
#>  310 :     std::vector<int> params_i;
#>  311 :     vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
#>  312 :              std::numeric_limits<double>::quiet_NaN());
#>  313 :     write_array_impl(base_rng, params_r, params_i, vars,
#>  314 :       emit_transformed_parameters, emit_generated_quantities, pstream);
#>  315 :   }
#>  316 :   template <typename RNG> inline void
#>  317 :   write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
#>  318 :               params_i, std::vector<double>& vars, bool
#>  319 :               emit_transformed_parameters = true, bool
#>  320 :               emit_generated_quantities = true, std::ostream*
#>  321 :               pstream = nullptr) const {
#>  322 :     const size_t num_params__ = 1;
#>  323 :     const size_t num_transformed = emit_transformed_parameters * (0U);
#>  324 :     const size_t num_gen_quantities = emit_generated_quantities * (0U);
#>  325 :     const size_t num_to_write = num_params__ + num_transformed +
#>  326 :       num_gen_quantities;
#>  327 :     vars = std::vector<double>(num_to_write,
#>  328 :              std::numeric_limits<double>::quiet_NaN());
#>  329 :     write_array_impl(base_rng, params_r, params_i, vars,
#>  330 :       emit_transformed_parameters, emit_generated_quantities, pstream);
#>  331 :   }
#>  332 :   template <bool propto__, bool jacobian__, typename T_> inline T_
#>  333 :   log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
#>  334 :     Eigen::Matrix<int,-1,1> params_i;
#>  335 :     return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
#>  336 :   }
#>  337 :   template <bool propto__, bool jacobian__, typename T_> inline T_
#>  338 :   log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
#>  339 :            std::ostream* pstream = nullptr) const {
#>  340 :     return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
#>  341 :   }
#>  342 :   inline void
#>  343 :   transform_inits(const stan::io::var_context& context,
#>  344 :                   Eigen::Matrix<double,-1,1>& params_r, std::ostream*
#>  345 :                   pstream = nullptr) const final {
#>  346 :     std::vector<double> params_r_vec(params_r.size());
#>  347 :     std::vector<int> params_i;
#>  348 :     transform_inits(context, params_i, params_r_vec, pstream);
#>  349 :     params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
#>  350 :                  params_r_vec.size());
#>  351 :   }
#>  352 :   inline void
#>  353 :   transform_inits(const stan::io::var_context& context, std::vector<int>&
#>  354 :                   params_i, std::vector<double>& vars, std::ostream*
#>  355 :                   pstream__ = nullptr) const {
#>  356 :     vars.resize(num_params_r__);
#>  357 :     transform_inits_impl(context, vars, pstream__);
#>  358 :   }
#>  359 :   inline void
#>  360 :   unconstrain_array(const std::vector<double>& params_constrained,
#>  361 :                     std::vector<double>& params_unconstrained, std::ostream*
#>  362 :                     pstream = nullptr) const {
#>  363 :     const std::vector<int> params_i;
#>  364 :     params_unconstrained = std::vector<double>(num_params_r__,
#>  365 :                              std::numeric_limits<double>::quiet_NaN());
#>  366 :     unconstrain_array_impl(params_constrained, params_i,
#>  367 :       params_unconstrained, pstream);
#>  368 :   }
#>  369 :   inline void
#>  370 :   unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
#>  371 :                     Eigen::Matrix<double,-1,1>& params_unconstrained,
#>  372 :                     std::ostream* pstream = nullptr) const {
#>  373 :     const std::vector<int> params_i;
#>  374 :     params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
#>  375 :                              std::numeric_limits<double>::quiet_NaN());
#>  376 :     unconstrain_array_impl(params_constrained, params_i,
#>  377 :       params_unconstrained, pstream);
#>  378 :   }
#>  379 : };
#>  380 : }
#>  381 : using stan_model = model1c8428b78e6e_example_namespace::model1c8428b78e6e_example;
#>  382 : #ifndef USING_R
#>  383 : // Boilerplate
#>  384 : stan::model::model_base&
#>  385 : new_model(stan::io::var_context& data_context, unsigned int seed,
#>  386 :           std::ostream* msg_stream) {
#>  387 :   stan_model* m = new stan_model(data_context, seed, msg_stream);
#>  388 :   return *m;
#>  389 : }
#>  390 : stan::math::profile_map& get_stan_profile_data() {
#>  391 :   return model1c8428b78e6e_example_namespace::profiles__;
#>  392 : }
#>  393 : #endif
#>  394 : #endif
#>  395 : 
#>  396 : RCPP_MODULE(stan_fit4model1c8428b78e6e_example_mod) {
#>  397 :   class_<rstan::stan_fit<stan_model, boost::random::mixmax> >(
#>  398 :       "stan_fit4model1c8428b78e6e_example")
#>  399 : 
#>  400 :       .constructor<SEXP, SEXP, SEXP>()
#>  401 : 
#>  402 :       .method(
#>  403 :           "call_sampler",
#>  404 :           &rstan::stan_fit<stan_model, boost::random::mixmax>::call_sampler)
#>  405 :       .method(
#>  406 :           "param_names",
#>  407 :           &rstan::stan_fit<stan_model, boost::random::mixmax>::param_names)
#>  408 :       .method("param_names_oi",
#>  409 :               &rstan::stan_fit<stan_model,
#>  410 :                                boost::random::mixmax>::param_names_oi)
#>  411 :       .method("param_fnames_oi",
#>  412 :               &rstan::stan_fit<stan_model,
#>  413 :                                boost::random::mixmax>::param_fnames_oi)
#>  414 :       .method(
#>  415 :           "param_dims",
#>  416 :           &rstan::stan_fit<stan_model, boost::random::mixmax>::param_dims)
#>  417 :       .method("param_dims_oi",
#>  418 :               &rstan::stan_fit<stan_model,
#>  419 :                                boost::random::mixmax>::param_dims_oi)
#>  420 :       .method("update_param_oi",
#>  421 :               &rstan::stan_fit<stan_model,
#>  422 :                                boost::random::mixmax>::update_param_oi)
#>  423 :       .method("param_oi_tidx",
#>  424 :               &rstan::stan_fit<stan_model,
#>  425 :                                boost::random::mixmax>::param_oi_tidx)
#>  426 :       .method("grad_log_prob",
#>  427 :               &rstan::stan_fit<stan_model,
#>  428 :                                boost::random::mixmax>::grad_log_prob)
#>  429 :       .method("log_prob",
#>  430 :               &rstan::stan_fit<stan_model, boost::random::mixmax>::log_prob)
#>  431 :       .method("unconstrain_pars",
#>  432 :               &rstan::stan_fit<stan_model,
#>  433 :                                boost::random::mixmax>::unconstrain_pars)
#>  434 :       .method("constrain_pars",
#>  435 :               &rstan::stan_fit<stan_model,
#>  436 :                                boost::random::mixmax>::constrain_pars)
#>  437 :       .method(
#>  438 :           "num_pars_unconstrained",
#>  439 :           &rstan::stan_fit<stan_model,
#>  440 :                            boost::random::mixmax>::num_pars_unconstrained)
#>  441 :       .method(
#>  442 :           "unconstrained_param_names",
#>  443 :           &rstan::stan_fit<
#>  444 :               stan_model, boost::random::mixmax>::unconstrained_param_names)
#>  445 :       .method(
#>  446 :           "constrained_param_names",
#>  447 :           &rstan::stan_fit<stan_model,
#>  448 :                            boost::random::mixmax>::constrained_param_names)
#>  449 :       .method("standalone_gqs",
#>  450 :               &rstan::stan_fit<stan_model,
#>  451 :                                boost::random::mixmax>::standalone_gqs);
#>  452 : }
#>  453 : 
#>  454 : 
#>  455 : // declarations
#>  456 : extern "C" {
#>  457 : SEXP file1c8412017029( ) ;
#>  458 : }
#>  459 : 
#>  460 : // definition
#>  461 : SEXP file1c8412017029() {
#>  462 :  return Rcpp::wrap("example");
#>  463 : }
#> 
#> CHECKING DATA AND PREPROCESSING FOR MODEL 'example' NOW.
#> 
#> COMPILING MODEL 'example' NOW.
#> 
#> STARTING SAMPLER FOR MODEL 'example' NOW.
#> 
#> SAMPLING FOR MODEL 'example' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 3e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.03 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2012 [  0%]  (Warmup)
#> Chain 1: Iteration:  201 / 2012 [  9%]  (Warmup)
#> Chain 1: Iteration:  402 / 2012 [ 19%]  (Warmup)
#> Chain 1: Iteration:  603 / 2012 [ 29%]  (Warmup)
#> Chain 1: Iteration:  804 / 2012 [ 39%]  (Warmup)
#> Chain 1: Iteration: 1005 / 2012 [ 49%]  (Warmup)
#> Chain 1: Iteration: 1007 / 2012 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1207 / 2012 [ 59%]  (Sampling)
#> Chain 1: Iteration: 1408 / 2012 [ 69%]  (Sampling)
#> Chain 1: Iteration: 1609 / 2012 [ 79%]  (Sampling)
#> Chain 1: Iteration: 1810 / 2012 [ 89%]  (Sampling)
#> Chain 1: Iteration: 2011 / 2012 [ 99%]  (Sampling)
#> Chain 1: Iteration: 2012 / 2012 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.007 seconds (Warm-up)
#> Chain 1:                0.007 seconds (Sampling)
#> Chain 1:                0.014 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'example' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.02 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2012 [  0%]  (Warmup)
#> Chain 2: Iteration:  201 / 2012 [  9%]  (Warmup)
#> Chain 2: Iteration:  402 / 2012 [ 19%]  (Warmup)
#> Chain 2: Iteration:  603 / 2012 [ 29%]  (Warmup)
#> Chain 2: Iteration:  804 / 2012 [ 39%]  (Warmup)
#> Chain 2: Iteration: 1005 / 2012 [ 49%]  (Warmup)
#> Chain 2: Iteration: 1007 / 2012 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1207 / 2012 [ 59%]  (Sampling)
#> Chain 2: Iteration: 1408 / 2012 [ 69%]  (Sampling)
#> Chain 2: Iteration: 1609 / 2012 [ 79%]  (Sampling)
#> Chain 2: Iteration: 1810 / 2012 [ 89%]  (Sampling)
#> Chain 2: Iteration: 2011 / 2012 [ 99%]  (Sampling)
#> Chain 2: Iteration: 2012 / 2012 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.007 seconds (Warm-up)
#> Chain 2:                0.007 seconds (Sampling)
#> Chain 2:                0.014 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'example' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 2e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.02 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2012 [  0%]  (Warmup)
#> Chain 3: Iteration:  201 / 2012 [  9%]  (Warmup)
#> Chain 3: Iteration:  402 / 2012 [ 19%]  (Warmup)
#> Chain 3: Iteration:  603 / 2012 [ 29%]  (Warmup)
#> Chain 3: Iteration:  804 / 2012 [ 39%]  (Warmup)
#> Chain 3: Iteration: 1005 / 2012 [ 49%]  (Warmup)
#> Chain 3: Iteration: 1007 / 2012 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1207 / 2012 [ 59%]  (Sampling)
#> Chain 3: Iteration: 1408 / 2012 [ 69%]  (Sampling)
#> Chain 3: Iteration: 1609 / 2012 [ 79%]  (Sampling)
#> Chain 3: Iteration: 1810 / 2012 [ 89%]  (Sampling)
#> Chain 3: Iteration: 2011 / 2012 [ 99%]  (Sampling)
#> Chain 3: Iteration: 2012 / 2012 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.007 seconds (Warm-up)
#> Chain 3:                0.007 seconds (Sampling)
#> Chain 3:                0.014 seconds (Total)
#> Chain 3: 
print(fit)
#> Inference for Stan model: example.
#> 3 chains, each with iter=2012; warmup=1006; thin=1; 
#> post-warmup draws per chain=1006, total post-warmup draws=3018.
#> 
#>        mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
#> mu    -0.06    0.01 0.22  -0.50  -0.21  -0.06   0.09   0.38  1094    1
#> lp__ -35.26    0.02 0.71 -37.35 -35.41 -34.98 -34.81 -34.76  1086    1
#> 
#> Samples were drawn using NUTS(diag_e) at Fri Jul 24 23:30:10 2026.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).

# extract samples
e <- extract(fit, permuted = FALSE) # return a list of arrays
str(e)
#>  num [1:1006, 1:3, 1:2] 0.1563 -0.0556 -0.0701 -0.3364 -0.1217 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ iterations: NULL
#>   ..$ chains    : chr [1:3] "chain:1" "chain:2" "chain:3"
#>   ..$ parameters: chr [1:2] "mu" "lp__"

arr <- as.array(fit) # return an array
str(arr)
#>  num [1:1006, 1:3, 1:2] 0.1563 -0.0556 -0.0701 -0.3364 -0.1217 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ iterations: NULL
#>   ..$ chains    : chr [1:3] "chain:1" "chain:2" "chain:3"
#>   ..$ parameters: chr [1:2] "mu" "lp__"

mat <- as.matrix(fit) # return a matrix
str(mat)
#>  num [1:3018, 1:2] 0.1563 -0.0556 -0.0701 -0.3364 -0.1217 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ iterations: NULL
#>   ..$ parameters: chr [1:2] "mu" "lp__"
# }
```
