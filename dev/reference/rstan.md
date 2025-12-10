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

|             |                                                                                             |                                                                                        |                                                           |
|-------------|---------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------|-----------------------------------------------------------|
| **Package** | **Description**                                                                             | **Doc**                                                                                | **Website**                                               |
| bayesplot   | ggplot-based plotting of parameter estimates, diagnostics, and posterior predictive checks. | [bayesplot-package](https://mc-stan.org/bayesplot/reference/bayesplot-package.html)    | [mc-stan.org/bayesplot](https://mc-stan.org/bayesplot/)   |
| shinystan   | Interactive GUI for exploring MCMC output.                                                  | [shinystan-package](https://mc-stan.org/shinystan/reference/shinystan-package.html)    | [mc-stan.org/shinystan](https://mc-stan.org/shinystan/)   |
| loo         | Out-of-sample predictive performance estimates and model comparison.                        | [loo-package](https://mc-stan.org/loo/reference/loo-package.html)                      | [mc-stan.org/loo](https://mc-stan.org/loo/)               |
| rstanarm    | R formula interface for applied regression modeling.                                        | rstanarm-package                                                                       | [mc-stan.org/rstanarm](https://mc-stan.org/rstanarm/)     |
| rstantools  | Tools for developers of R packages interfacing with Stan.                                   | [rstantools-package](https://mc-stan.org/rstantools/reference/rstantools-package.html) | [mc-stan.org/rstantools](https://mc-stan.org/rstantools/) |

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
#> OS: x86_64, linux-gnu; rstan: 2.36.0.9000; Rcpp: 1.1.0; inline: 0.3.21 
#>  >> setting environment variables: 
#> PKG_LIBS =  '/home/runner/work/_temp/Library/rstan/lib//libStanServices.a' -L'/home/runner/work/_temp/Library/StanHeaders/lib/' -lStanHeaders -L'/home/runner/work/_temp/Library/RcppParallel/lib/' -ltbb 
#> PKG_CPPFLAGS =   -I"/home/runner/work/_temp/Library/Rcpp/include/"  -I"/home/runner/work/_temp/Library/RcppEigen/include/"  -I"/home/runner/work/_temp/Library/RcppEigen/include/unsupported"  -I"/home/runner/work/_temp/Library/BH/include" -I"/home/runner/work/_temp/Library/StanHeaders/include/src/"  -I"/home/runner/work/_temp/Library/StanHeaders/include/"  -I"/home/runner/work/_temp/Library/RcppParallel/include/"  -I"/home/runner/work/_temp/Library/rstan/include" -DEIGEN_NO_DEBUG  -DBOOST_DISABLE_ASSERTS  -DBOOST_PENDING_INTEGER_LOG2_HPP  -DSTAN_THREADS  -DUSE_STANC3 -DSTRICT_R_HEADERS  -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION  -D_HAS_AUTO_PTR_ETC=0  -include '/home/runner/work/_temp/Library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp'  -D_REENTRANT -DRCPP_PARALLEL_USE_TBB=1
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
#>   16 : // Code generated by stanc v2.37.0-172-g8d78669
#>   17 : #include <stan/model/model_header.hpp>
#>   18 : namespace model1bf852f26221_example_namespace {
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
#>   30 : class model1bf852f26221_example final : public model_base_crtp<model1bf852f26221_example> {
#>   31 : private:
#>   32 :   int N;
#>   33 :   std::vector<double> y;
#>   34 : public:
#>   35 :   ~model1bf852f26221_example() {}
#>   36 :   model1bf852f26221_example(stan::io::var_context& context__, unsigned int
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
#>   47 :       "model1bf852f26221_example_namespace::model1bf852f26221_example";
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
#>   76 :     return "model1bf852f26221_example";
#>   77 :   }
#>   78 :   inline std::vector<std::string> model_compile_info() const noexcept {
#>   79 :     return std::vector<std::string>{"stanc_version = stanc3 v2.37.0-172-g8d78669",
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
#>  102 :       "model1bf852f26221_example_namespace::log_prob";
#>  103 :     // suppress unused var warning
#>  104 :     (void) function__;
#>  105 :     try {
#>  106 :       current_statement__ = 1;
#>  107 :       auto mu = in__.template read<local_scalar_t__>();
#>  108 :       {
#>  109 :         current_statement__ = 2;
#>  110 :         lp_accum__.add(stan::math::normal_lpdf<false>(mu, 0, 10));
#>  111 :         current_statement__ = 3;
#>  112 :         lp_accum__.add(stan::math::normal_lpdf<false>(y, mu, 1));
#>  113 :       }
#>  114 :     } catch (const std::exception& e) {
#>  115 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  116 :     }
#>  117 :     lp_accum__.add(lp__);
#>  118 :     return lp_accum__.sum();
#>  119 :   }
#>  120 :   // Reverse mode autodiff log prob
#>  121 :   template <bool propto__, bool jacobian__, typename VecR, typename VecI,
#>  122 :             stan::require_vector_like_t<VecR>* = nullptr,
#>  123 :             stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr,
#>  124 :             stan::require_st_var<VecR>* = nullptr>
#>  125 :   inline stan::scalar_type_t<VecR>
#>  126 :   log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
#>  127 :                 pstream__ = nullptr) const {
#>  128 :     using T__ = stan::scalar_type_t<VecR>;
#>  129 :     using local_scalar_t__ = T__;
#>  130 :     T__ lp__(0.0);
#>  131 :     stan::math::accumulator<T__> lp_accum__;
#>  132 :     stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
#>  133 :     int current_statement__ = 0;
#>  134 :     // suppress unused var warning
#>  135 :     (void) current_statement__;
#>  136 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>  137 :     // suppress unused var warning
#>  138 :     (void) DUMMY_VAR__;
#>  139 :     static constexpr const char* function__ =
#>  140 :       "model1bf852f26221_example_namespace::log_prob";
#>  141 :     // suppress unused var warning
#>  142 :     (void) function__;
#>  143 :     try {
#>  144 :       current_statement__ = 1;
#>  145 :       auto mu = in__.template read<local_scalar_t__>();
#>  146 :       {
#>  147 :         current_statement__ = 2;
#>  148 :         lp_accum__.add(stan::math::normal_lpdf<false>(mu, 0, 10));
#>  149 :         current_statement__ = 3;
#>  150 :         lp_accum__.add(stan::math::normal_lpdf<false>(y, mu, 1));
#>  151 :       }
#>  152 :     } catch (const std::exception& e) {
#>  153 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  154 :     }
#>  155 :     lp_accum__.add(lp__);
#>  156 :     return lp_accum__.sum();
#>  157 :   }
#>  158 :   template <typename RNG, typename VecR, typename VecI, typename VecVar,
#>  159 :             stan::require_vector_like_vt<std::is_floating_point,
#>  160 :             VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
#>  161 :             VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
#>  162 :             VecVar>* = nullptr>
#>  163 :   inline void
#>  164 :   write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
#>  165 :                    VecVar& vars__, const bool
#>  166 :                    emit_transformed_parameters__ = true, const bool
#>  167 :                    emit_generated_quantities__ = true, std::ostream*
#>  168 :                    pstream__ = nullptr) const {
#>  169 :     using local_scalar_t__ = double;
#>  170 :     stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
#>  171 :     stan::io::serializer<local_scalar_t__> out__(vars__);
#>  172 :     static constexpr bool propto__ = true;
#>  173 :     // suppress unused var warning
#>  174 :     (void) propto__;
#>  175 :     double lp__ = 0.0;
#>  176 :     // suppress unused var warning
#>  177 :     (void) lp__;
#>  178 :     int current_statement__ = 0;
#>  179 :     // suppress unused var warning
#>  180 :     (void) current_statement__;
#>  181 :     stan::math::accumulator<double> lp_accum__;
#>  182 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>  183 :     // suppress unused var warning
#>  184 :     (void) DUMMY_VAR__;
#>  185 :     constexpr bool jacobian__ = false;
#>  186 :     // suppress unused var warning
#>  187 :     (void) jacobian__;
#>  188 :     static constexpr const char* function__ =
#>  189 :       "model1bf852f26221_example_namespace::write_array";
#>  190 :     // suppress unused var warning
#>  191 :     (void) function__;
#>  192 :     try {
#>  193 :       current_statement__ = 1;
#>  194 :       auto mu = in__.template read<local_scalar_t__>();
#>  195 :       out__.write(mu);
#>  196 :       if (stan::math::logical_negation(
#>  197 :             (stan::math::primitive_value(emit_transformed_parameters__) ||
#>  198 :             stan::math::primitive_value(emit_generated_quantities__)))) {
#>  199 :         return ;
#>  200 :       }
#>  201 :       if (stan::math::logical_negation(emit_generated_quantities__)) {
#>  202 :         return ;
#>  203 :       }
#>  204 :     } catch (const std::exception& e) {
#>  205 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  206 :     }
#>  207 :   }
#>  208 :   template <typename VecVar, typename VecI,
#>  209 :             stan::require_vector_t<VecVar>* = nullptr,
#>  210 :             stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
#>  211 :   inline void
#>  212 :   unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
#>  213 :                          VecVar& vars__, std::ostream* pstream__ = nullptr) const {
#>  214 :     using local_scalar_t__ = double;
#>  215 :     stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
#>  216 :     stan::io::serializer<local_scalar_t__> out__(vars__);
#>  217 :     int current_statement__ = 0;
#>  218 :     // suppress unused var warning
#>  219 :     (void) current_statement__;
#>  220 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>  221 :     // suppress unused var warning
#>  222 :     (void) DUMMY_VAR__;
#>  223 :     try {
#>  224 :       local_scalar_t__ mu = DUMMY_VAR__;
#>  225 :       current_statement__ = 1;
#>  226 :       mu = in__.read<local_scalar_t__>();
#>  227 :       out__.write(mu);
#>  228 :     } catch (const std::exception& e) {
#>  229 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  230 :     }
#>  231 :   }
#>  232 :   template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
#>  233 :   inline void
#>  234 :   transform_inits_impl(const stan::io::var_context& context__, VecVar&
#>  235 :                        vars__, std::ostream* pstream__ = nullptr) const {
#>  236 :     using local_scalar_t__ = double;
#>  237 :     stan::io::serializer<local_scalar_t__> out__(vars__);
#>  238 :     int current_statement__ = 0;
#>  239 :     // suppress unused var warning
#>  240 :     (void) current_statement__;
#>  241 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>  242 :     // suppress unused var warning
#>  243 :     (void) DUMMY_VAR__;
#>  244 :     try {
#>  245 :       current_statement__ = 1;
#>  246 :       context__.validate_dims("parameter initialization", "mu", "double",
#>  247 :         std::vector<size_t>{});
#>  248 :       local_scalar_t__ mu = DUMMY_VAR__;
#>  249 :       current_statement__ = 1;
#>  250 :       mu = context__.vals_r("mu")[(1 - 1)];
#>  251 :       out__.write(mu);
#>  252 :     } catch (const std::exception& e) {
#>  253 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  254 :     }
#>  255 :   }
#>  256 :   inline void
#>  257 :   get_param_names(std::vector<std::string>& names__, const bool
#>  258 :                   emit_transformed_parameters__ = true, const bool
#>  259 :                   emit_generated_quantities__ = true) const {
#>  260 :     names__ = std::vector<std::string>{"mu"};
#>  261 :     if (emit_transformed_parameters__) {}
#>  262 :     if (emit_generated_quantities__) {}
#>  263 :   }
#>  264 :   inline void
#>  265 :   get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
#>  266 :            emit_transformed_parameters__ = true, const bool
#>  267 :            emit_generated_quantities__ = true) const {
#>  268 :     dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{}};
#>  269 :     if (emit_transformed_parameters__) {}
#>  270 :     if (emit_generated_quantities__) {}
#>  271 :   }
#>  272 :   inline void
#>  273 :   constrained_param_names(std::vector<std::string>& param_names__, bool
#>  274 :                           emit_transformed_parameters__ = true, bool
#>  275 :                           emit_generated_quantities__ = true) const final {
#>  276 :     param_names__.emplace_back(std::string() + "mu");
#>  277 :     if (emit_transformed_parameters__) {}
#>  278 :     if (emit_generated_quantities__) {}
#>  279 :   }
#>  280 :   inline void
#>  281 :   unconstrained_param_names(std::vector<std::string>& param_names__, bool
#>  282 :                             emit_transformed_parameters__ = true, bool
#>  283 :                             emit_generated_quantities__ = true) const final {
#>  284 :     param_names__.emplace_back(std::string() + "mu");
#>  285 :     if (emit_transformed_parameters__) {}
#>  286 :     if (emit_generated_quantities__) {}
#>  287 :   }
#>  288 :   inline std::string get_constrained_sizedtypes() const {
#>  289 :     return std::string("[{\"name\":\"mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]");
#>  290 :   }
#>  291 :   inline std::string get_unconstrained_sizedtypes() const {
#>  292 :     return std::string("[{\"name\":\"mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]");
#>  293 :   }
#>  294 :   // Begin method overload boilerplate
#>  295 :   template <typename RNG> inline void
#>  296 :   write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
#>  297 :               Eigen::Matrix<double,-1,1>& vars, const bool
#>  298 :               emit_transformed_parameters = true, const bool
#>  299 :               emit_generated_quantities = true, std::ostream*
#>  300 :               pstream = nullptr) const {
#>  301 :     const size_t num_params__ = 1;
#>  302 :     const size_t num_transformed = emit_transformed_parameters * (0U);
#>  303 :     const size_t num_gen_quantities = emit_generated_quantities * (0U);
#>  304 :     const size_t num_to_write = num_params__ + num_transformed +
#>  305 :       num_gen_quantities;
#>  306 :     std::vector<int> params_i;
#>  307 :     vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
#>  308 :              std::numeric_limits<double>::quiet_NaN());
#>  309 :     write_array_impl(base_rng, params_r, params_i, vars,
#>  310 :       emit_transformed_parameters, emit_generated_quantities, pstream);
#>  311 :   }
#>  312 :   template <typename RNG> inline void
#>  313 :   write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
#>  314 :               params_i, std::vector<double>& vars, bool
#>  315 :               emit_transformed_parameters = true, bool
#>  316 :               emit_generated_quantities = true, std::ostream*
#>  317 :               pstream = nullptr) const {
#>  318 :     const size_t num_params__ = 1;
#>  319 :     const size_t num_transformed = emit_transformed_parameters * (0U);
#>  320 :     const size_t num_gen_quantities = emit_generated_quantities * (0U);
#>  321 :     const size_t num_to_write = num_params__ + num_transformed +
#>  322 :       num_gen_quantities;
#>  323 :     vars = std::vector<double>(num_to_write,
#>  324 :              std::numeric_limits<double>::quiet_NaN());
#>  325 :     write_array_impl(base_rng, params_r, params_i, vars,
#>  326 :       emit_transformed_parameters, emit_generated_quantities, pstream);
#>  327 :   }
#>  328 :   template <bool propto__, bool jacobian__, typename T_> inline T_
#>  329 :   log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
#>  330 :     Eigen::Matrix<int,-1,1> params_i;
#>  331 :     return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
#>  332 :   }
#>  333 :   template <bool propto__, bool jacobian__, typename T_> inline T_
#>  334 :   log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
#>  335 :            std::ostream* pstream = nullptr) const {
#>  336 :     return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
#>  337 :   }
#>  338 :   inline void
#>  339 :   transform_inits(const stan::io::var_context& context,
#>  340 :                   Eigen::Matrix<double,-1,1>& params_r, std::ostream*
#>  341 :                   pstream = nullptr) const final {
#>  342 :     std::vector<double> params_r_vec(params_r.size());
#>  343 :     std::vector<int> params_i;
#>  344 :     transform_inits(context, params_i, params_r_vec, pstream);
#>  345 :     params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
#>  346 :                  params_r_vec.size());
#>  347 :   }
#>  348 :   inline void
#>  349 :   transform_inits(const stan::io::var_context& context, std::vector<int>&
#>  350 :                   params_i, std::vector<double>& vars, std::ostream*
#>  351 :                   pstream__ = nullptr) const {
#>  352 :     vars.resize(num_params_r__);
#>  353 :     transform_inits_impl(context, vars, pstream__);
#>  354 :   }
#>  355 :   inline void
#>  356 :   unconstrain_array(const std::vector<double>& params_constrained,
#>  357 :                     std::vector<double>& params_unconstrained, std::ostream*
#>  358 :                     pstream = nullptr) const {
#>  359 :     const std::vector<int> params_i;
#>  360 :     params_unconstrained = std::vector<double>(num_params_r__,
#>  361 :                              std::numeric_limits<double>::quiet_NaN());
#>  362 :     unconstrain_array_impl(params_constrained, params_i,
#>  363 :       params_unconstrained, pstream);
#>  364 :   }
#>  365 :   inline void
#>  366 :   unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
#>  367 :                     Eigen::Matrix<double,-1,1>& params_unconstrained,
#>  368 :                     std::ostream* pstream = nullptr) const {
#>  369 :     const std::vector<int> params_i;
#>  370 :     params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
#>  371 :                              std::numeric_limits<double>::quiet_NaN());
#>  372 :     unconstrain_array_impl(params_constrained, params_i,
#>  373 :       params_unconstrained, pstream);
#>  374 :   }
#>  375 : };
#>  376 : }
#>  377 : using stan_model = model1bf852f26221_example_namespace::model1bf852f26221_example;
#>  378 : #ifndef USING_R
#>  379 : // Boilerplate
#>  380 : stan::model::model_base&
#>  381 : new_model(stan::io::var_context& data_context, unsigned int seed,
#>  382 :           std::ostream* msg_stream) {
#>  383 :   stan_model* m = new stan_model(data_context, seed, msg_stream);
#>  384 :   return *m;
#>  385 : }
#>  386 : stan::math::profile_map& get_stan_profile_data() {
#>  387 :   return model1bf852f26221_example_namespace::profiles__;
#>  388 : }
#>  389 : #endif
#>  390 : #endif
#>  391 : 
#>  392 : RCPP_MODULE(stan_fit4model1bf852f26221_example_mod) {
#>  393 :   class_<rstan::stan_fit<stan_model, boost::random::mixmax> >(
#>  394 :       "stan_fit4model1bf852f26221_example")
#>  395 : 
#>  396 :       .constructor<SEXP, SEXP, SEXP>()
#>  397 : 
#>  398 :       .method(
#>  399 :           "call_sampler",
#>  400 :           &rstan::stan_fit<stan_model, boost::random::mixmax>::call_sampler)
#>  401 :       .method(
#>  402 :           "param_names",
#>  403 :           &rstan::stan_fit<stan_model, boost::random::mixmax>::param_names)
#>  404 :       .method("param_names_oi",
#>  405 :               &rstan::stan_fit<stan_model,
#>  406 :                                boost::random::mixmax>::param_names_oi)
#>  407 :       .method("param_fnames_oi",
#>  408 :               &rstan::stan_fit<stan_model,
#>  409 :                                boost::random::mixmax>::param_fnames_oi)
#>  410 :       .method(
#>  411 :           "param_dims",
#>  412 :           &rstan::stan_fit<stan_model, boost::random::mixmax>::param_dims)
#>  413 :       .method("param_dims_oi",
#>  414 :               &rstan::stan_fit<stan_model,
#>  415 :                                boost::random::mixmax>::param_dims_oi)
#>  416 :       .method("update_param_oi",
#>  417 :               &rstan::stan_fit<stan_model,
#>  418 :                                boost::random::mixmax>::update_param_oi)
#>  419 :       .method("param_oi_tidx",
#>  420 :               &rstan::stan_fit<stan_model,
#>  421 :                                boost::random::mixmax>::param_oi_tidx)
#>  422 :       .method("grad_log_prob",
#>  423 :               &rstan::stan_fit<stan_model,
#>  424 :                                boost::random::mixmax>::grad_log_prob)
#>  425 :       .method("log_prob",
#>  426 :               &rstan::stan_fit<stan_model, boost::random::mixmax>::log_prob)
#>  427 :       .method("unconstrain_pars",
#>  428 :               &rstan::stan_fit<stan_model,
#>  429 :                                boost::random::mixmax>::unconstrain_pars)
#>  430 :       .method("constrain_pars",
#>  431 :               &rstan::stan_fit<stan_model,
#>  432 :                                boost::random::mixmax>::constrain_pars)
#>  433 :       .method(
#>  434 :           "num_pars_unconstrained",
#>  435 :           &rstan::stan_fit<stan_model,
#>  436 :                            boost::random::mixmax>::num_pars_unconstrained)
#>  437 :       .method(
#>  438 :           "unconstrained_param_names",
#>  439 :           &rstan::stan_fit<
#>  440 :               stan_model, boost::random::mixmax>::unconstrained_param_names)
#>  441 :       .method(
#>  442 :           "constrained_param_names",
#>  443 :           &rstan::stan_fit<stan_model,
#>  444 :                            boost::random::mixmax>::constrained_param_names)
#>  445 :       .method("standalone_gqs",
#>  446 :               &rstan::stan_fit<stan_model,
#>  447 :                                boost::random::mixmax>::standalone_gqs);
#>  448 : }
#>  449 : 
#>  450 : 
#>  451 : // declarations
#>  452 : extern "C" {
#>  453 : SEXP file1bf84a9e46df( ) ;
#>  454 : }
#>  455 : 
#>  456 : // definition
#>  457 : SEXP file1bf84a9e46df() {
#>  458 :  return Rcpp::wrap("example");
#>  459 : }
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
#> Chain 2: Gradient evaluation took 1e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.01 seconds.
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
#> Chain 3: Gradient evaluation took 1e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.01 seconds.
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
#> mu    -0.22    0.01 0.23  -0.68  -0.37  -0.22  -0.07   0.21  1204    1
#> lp__ -30.62    0.02 0.75 -32.67 -30.76 -30.33 -30.15 -30.10  1058    1
#> 
#> Samples were drawn using NUTS(diag_e) at Wed Dec 10 00:46:30 2025.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).

# extract samples
e <- extract(fit, permuted = FALSE) # return a list of arrays
str(e)
#>  num [1:1006, 1:3, 1:2] -0.211 -0.272 -0.152 -0.24 0.15 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ iterations: NULL
#>   ..$ chains    : chr [1:3] "chain:1" "chain:2" "chain:3"
#>   ..$ parameters: chr [1:2] "mu" "lp__"

arr <- as.array(fit) # return an array
str(arr)
#>  num [1:1006, 1:3, 1:2] -0.211 -0.272 -0.152 -0.24 0.15 ...
#>  - attr(*, "dimnames")=List of 3
#>   ..$ iterations: NULL
#>   ..$ chains    : chr [1:3] "chain:1" "chain:2" "chain:3"
#>   ..$ parameters: chr [1:2] "mu" "lp__"

mat <- as.matrix(fit) # return a matrix
str(mat)
#>  num [1:3018, 1:2] -0.211 -0.272 -0.152 -0.24 0.15 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ iterations: NULL
#>   ..$ parameters: chr [1:2] "mu" "lp__"
# }
```
