# Construct a Stan model

Construct an instance of S4 class `stanmodel` from a model specified in
Stan's modeling language. A `stanmodel` object can then be used to draw
samples from the model. The Stan program (the model expressed in the
Stan modeling language) is first translated to C++ code and then the C++
code for the model plus other auxiliary code is compiled into a dynamic
shared object (DSO) and then loaded. The loaded DSO for the model can be
executed to draw samples, allowing inference to be performed for the
model and data.

## Usage

``` r
stan_model(
    file, model_name = "anon_model",
    model_code = "", stanc_ret = NULL,
    boost_lib = NULL, eigen_lib = NULL,
    save_dso = TRUE, verbose = FALSE,
    auto_write = rstan_options("auto_write"),
    obfuscate_model_name = TRUE,
    allow_undefined = isTRUE(getOption("stanc.allow_undefined", FALSE)),
    allow_optimizations = isTRUE(getOption("stanc.allow_optimizations", FALSE)),
    standalone_functions = isTRUE(getOption("stanc.standalone_functions", FALSE)),
    use_opencl = isTRUE(getOption("stanc.use_opencl", FALSE)),
    warn_pedantic = isTRUE(getOption("stanc.warn_pedantic", FALSE)),
    warn_uninitialized = isTRUE(getOption("stanc.warn_uninitialized", FALSE)),
    includes = NULL,
    isystem = c(if (!missing(file)) dirname(file), getwd()))
```

## Arguments

- file:

  A character string or a connection that R supports specifying the Stan
  model specification in Stan's modeling language.

- model_name:

  A character string naming the model; defaults to `"anon_model"`.
  However, the model name will be derived from `file` or `model_code`
  (if `model_code` is the name of a character string object) if
  `model_name` is not specified.

- model_code:

  Either a character string containing the model specification or the
  name of a character string object in the workspace. This is an
  alternative to specifying the model via the `file` or `stanc_ret`
  arguments.

- stanc_ret:

  A named list returned from a previous call to the
  [`stanc`](https://mc-stan.org/rstan/dev/reference/stanc.md) function.
  The list can be used to specify the model instead of using the `file`
  or `model_code` arguments.

- boost_lib:

  The path to a version of the Boost C++ library to use instead of the
  one in the BH package.

- eigen_lib:

  The path to a version of the Eigen C++ library to use instead of the
  one in the RcppEigen package.

- save_dso:

  Logical, defaulting to `TRUE`, indicating whether the dynamic shared
  object (DSO) compiled from the C++ code for the model will be saved or
  not. If `TRUE`, we can draw samples from the same model in another R
  session using the saved DSO (i.e., without compiling the C++ code
  again).

- verbose:

  Logical, defaulting to `FALSE`, indicating whether to report
  additional intermediate output to the console, which might be helpful
  for debugging.

- auto_write:

  Logical, defaulting to the value of `rstan_options("auto_write")`,
  indicating whether to write the object to the hard disk using
  [`saveRDS`](https://rdrr.io/r/base/readRDS.html). Although this
  argument is `FALSE` by default, we recommend calling
  `rstan_options("auto_write" = TRUE)` in order to avoid unnecessary
  recompilations. If `file` is supplied and its
  [`dirname`](https://rdrr.io/r/base/basename.html) is writable, then
  the object will be written to that same directory, substituting a
  `.rds` extension for the `.stan` extension. Otherwise, the object will
  be written to the [`tempdir`](https://rdrr.io/r/base/tempfile.html).

- obfuscate_model_name:

  A logical scalar that is `TRUE` by default and passed to
  [`stanc`](https://mc-stan.org/rstan/dev/reference/stanc.md).

- allow_undefined:

  A logical scalar that is `FALSE` by default and passed to
  [`stanc`](https://mc-stan.org/rstan/dev/reference/stanc.md).

- allow_optimizations:

  A logical scalar that is `FALSE` by default and passed to
  [`stanc`](https://mc-stan.org/rstan/dev/reference/stanc.md).

- standalone_functions:

  A logical scalar that is `FALSE` by default and passed to
  [`stanc`](https://mc-stan.org/rstan/dev/reference/stanc.md).

- use_opencl:

  A logical scalar that is `FALSE` by default and passed to
  [`stanc`](https://mc-stan.org/rstan/dev/reference/stanc.md).

- warn_pedantic:

  A logical scalar that is `FALSE` by default and passed to
  [`stanc`](https://mc-stan.org/rstan/dev/reference/stanc.md).

- warn_uninitialized:

  A logical scalar that is `FALSE` by default and passed to
  [`stanc`](https://mc-stan.org/rstan/dev/reference/stanc.md).

- includes:

  If not `NULL` (the default), then a character vector of length one
  (possibly containing one or more `"\n"`) of the form
  `'#include "/full/path/to/my_header.hpp"'`, which will be inserted
  into the C++ code in the model's namespace and can be used to provide
  definitions of functions that are declared but not defined in `file`
  or `model_code` when `allow_undefined = TRUE`

- isystem:

  A character vector naming a path to look for file paths in `file` that
  are to be included within the Stan program named by `file`. See the
  Details section below.

## Details

If a previously compiled `stanmodel` exists on the hard drive, its
validity is checked and then returned without recompiling. The most
common form of invalidity seems to be Stan code that ends with a `}`
rather than a blank line, which causes the hash checker to think that
the current model is different than the one saved on the hard drive. To
avoid reading previously compiled `stanmodel`s from the hard drive,
supply the `stanc_ret` argument rather than the `file` or `model_code`
arguments.

There are three ways to specify the model's code for `stan_model`:

1.  parameter `model_code`: a character string containing the Stan model
    specification,

2.  parameter `file`: a file name (or a connection) from which to read
    the Stan model specification, or

3.  parameter `stanc_ret`: a list returned by `stanc` to be reused.

## Value

An instance of S4 class
[`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md)
that can be passed to the
[`sampling`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-sampling.md),
[`optimizing`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-optimizing.md),
and
[`vb`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-vb.md)
functions.

## References

The Stan Development Team *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org/>.

## See also

[`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md)
for details on the class.

[`sampling`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-sampling.md),
[`optimizing`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-optimizing.md),
and
[`vb`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-vb.md),
which take a `stanmodel` object as input, for estimating the model
parameters.

More details on Stan, including the full user's guide and reference
manual, can be found at <https://mc-stan.org/>.

## Examples

``` r
# \dontrun{
stancode <- 'data {real y_mean;} parameters {real y;} model {y ~ normal(y_mean,1);}'
mod <- stan_model(model_code = stancode, verbose = TRUE)
#> 
#> TRANSLATING MODEL '' FROM Stan CODE TO C++ CODE NOW.
#> OS: x86_64, linux-gnu; rstan: 2.36.0.9000; Rcpp: 1.1.1; inline: 0.3.21 
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
#>   16 : // Code generated by stanc v2.38.0-89-g5d63971
#>   17 : #include <stan/model/model_header.hpp>
#>   18 : namespace model1fb45d036b29__namespace {
#>   19 : using stan::model::model_base_crtp;
#>   20 : using namespace stan::math;
#>   21 : stan::math::profile_map profiles__;
#>   22 : static constexpr std::array<const char*, 4> locations_array__ =
#>   23 :   {" (found before start of program)",
#>   24 :   " (in 'anon_model', line 1, column 32 to column 39)",
#>   25 :   " (in 'anon_model', line 1, column 48 to column 69)",
#>   26 :   " (in 'anon_model', line 1, column 6 to column 18)"};
#>   27 : class model1fb45d036b29_ final : public model_base_crtp<model1fb45d036b29_> {
#>   28 : private:
#>   29 :   double y_mean;
#>   30 : public:
#>   31 :   ~model1fb45d036b29_() {}
#>   32 :   model1fb45d036b29_(stan::io::var_context& context__, unsigned int
#>   33 :                      random_seed__ = 0, std::ostream* pstream__ = nullptr)
#>   34 :       : model_base_crtp(0) {
#>   35 :     int current_statement__ = 0;
#>   36 :     // suppress unused var warning
#>   37 :     (void) current_statement__;
#>   38 :     using local_scalar_t__ = double;
#>   39 :     auto base_rng__ = stan::services::util::create_rng(random_seed__, 0);
#>   40 :     // suppress unused var warning
#>   41 :     (void) base_rng__;
#>   42 :     static constexpr const char* function__ =
#>   43 :       "model1fb45d036b29__namespace::model1fb45d036b29_";
#>   44 :     // suppress unused var warning
#>   45 :     (void) function__;
#>   46 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>   47 :     // suppress unused var warning
#>   48 :     (void) DUMMY_VAR__;
#>   49 :     try {
#>   50 :       current_statement__ = 3;
#>   51 :       context__.validate_dims("data initialization", "y_mean", "double",
#>   52 :         std::vector<size_t>{});
#>   53 :       y_mean = std::numeric_limits<double>::quiet_NaN();
#>   54 :       current_statement__ = 3;
#>   55 :       y_mean = context__.vals_r("y_mean")[(1 - 1)];
#>   56 :     } catch (const std::exception& e) {
#>   57 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>   58 :     }
#>   59 :     num_params_r__ = 1;
#>   60 :   }
#>   61 :   inline std::string model_name() const final {
#>   62 :     return "model1fb45d036b29_";
#>   63 :   }
#>   64 :   inline std::vector<std::string> model_compile_info() const noexcept {
#>   65 :     return std::vector<std::string>{"stanc_version = stanc3 v2.38.0-89-g5d63971",
#>   66 :              "stancflags = --"};
#>   67 :   }
#>   68 :   // Base log prob
#>   69 :   template <bool propto__, bool jacobian__, typename VecR, typename VecI,
#>   70 :             stan::require_vector_like_t<VecR>* = nullptr,
#>   71 :             stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr,
#>   72 :             stan::require_not_st_var<VecR>* = nullptr>
#>   73 :   inline stan::scalar_type_t<VecR>
#>   74 :   log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
#>   75 :                 pstream__ = nullptr) const {
#>   76 :     using T__ = stan::scalar_type_t<VecR>;
#>   77 :     using local_scalar_t__ = T__;
#>   78 :     T__ lp__(0.0);
#>   79 :     stan::math::accumulator<T__> lp_accum__;
#>   80 :     stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
#>   81 :     int current_statement__ = 0;
#>   82 :     // suppress unused var warning
#>   83 :     (void) current_statement__;
#>   84 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>   85 :     // suppress unused var warning
#>   86 :     (void) DUMMY_VAR__;
#>   87 :     static constexpr const char* function__ =
#>   88 :       "model1fb45d036b29__namespace::log_prob";
#>   89 :     // suppress unused var warning
#>   90 :     (void) function__;
#>   91 :     try {
#>   92 :       current_statement__ = 1;
#>   93 :       auto y = in__.template read<local_scalar_t__>();
#>   94 :       {
#>   95 :         current_statement__ = 2;
#>   96 :         lp_accum__.add(stan::math::normal_lpdf<propto__>(y, y_mean,
#>   97 :                          static_cast<double>(1)));
#>   98 :       }
#>   99 :     } catch (const std::exception& e) {
#>  100 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  101 :     }
#>  102 :     lp_accum__.add(lp__);
#>  103 :     return lp_accum__.sum();
#>  104 :   }
#>  105 :   // Reverse mode autodiff log prob
#>  106 :   template <bool propto__, bool jacobian__, typename VecR, typename VecI,
#>  107 :             stan::require_vector_like_t<VecR>* = nullptr,
#>  108 :             stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr,
#>  109 :             stan::require_st_var<VecR>* = nullptr>
#>  110 :   inline stan::scalar_type_t<VecR>
#>  111 :   log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
#>  112 :                 pstream__ = nullptr) const {
#>  113 :     using T__ = stan::scalar_type_t<VecR>;
#>  114 :     using local_scalar_t__ = T__;
#>  115 :     T__ lp__(0.0);
#>  116 :     stan::math::accumulator<T__> lp_accum__;
#>  117 :     stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
#>  118 :     int current_statement__ = 0;
#>  119 :     // suppress unused var warning
#>  120 :     (void) current_statement__;
#>  121 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>  122 :     // suppress unused var warning
#>  123 :     (void) DUMMY_VAR__;
#>  124 :     static constexpr const char* function__ =
#>  125 :       "model1fb45d036b29__namespace::log_prob";
#>  126 :     // suppress unused var warning
#>  127 :     (void) function__;
#>  128 :     try {
#>  129 :       current_statement__ = 1;
#>  130 :       auto y = in__.template read<local_scalar_t__>();
#>  131 :       {
#>  132 :         current_statement__ = 2;
#>  133 :         lp_accum__.add(stan::math::normal_lpdf<propto__>(y, y_mean,
#>  134 :                          static_cast<double>(1)));
#>  135 :       }
#>  136 :     } catch (const std::exception& e) {
#>  137 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  138 :     }
#>  139 :     lp_accum__.add(lp__);
#>  140 :     return lp_accum__.sum();
#>  141 :   }
#>  142 :   template <typename RNG, typename VecR, typename VecI, typename VecVar,
#>  143 :             stan::require_vector_like_vt<std::is_floating_point,
#>  144 :             VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
#>  145 :             VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
#>  146 :             VecVar>* = nullptr>
#>  147 :   inline void
#>  148 :   write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
#>  149 :                    VecVar& vars__, const bool
#>  150 :                    emit_transformed_parameters__ = true, const bool
#>  151 :                    emit_generated_quantities__ = true, std::ostream*
#>  152 :                    pstream__ = nullptr) const {
#>  153 :     using local_scalar_t__ = double;
#>  154 :     stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
#>  155 :     stan::io::serializer<local_scalar_t__> out__(vars__);
#>  156 :     static constexpr bool propto__ = true;
#>  157 :     // suppress unused var warning
#>  158 :     (void) propto__;
#>  159 :     double lp__ = 0.0;
#>  160 :     // suppress unused var warning
#>  161 :     (void) lp__;
#>  162 :     int current_statement__ = 0;
#>  163 :     // suppress unused var warning
#>  164 :     (void) current_statement__;
#>  165 :     stan::math::accumulator<double> lp_accum__;
#>  166 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>  167 :     // suppress unused var warning
#>  168 :     (void) DUMMY_VAR__;
#>  169 :     constexpr bool jacobian__ = false;
#>  170 :     // suppress unused var warning
#>  171 :     (void) jacobian__;
#>  172 :     static constexpr const char* function__ =
#>  173 :       "model1fb45d036b29__namespace::write_array";
#>  174 :     // suppress unused var warning
#>  175 :     (void) function__;
#>  176 :     try {
#>  177 :       current_statement__ = 1;
#>  178 :       auto y = in__.template read<local_scalar_t__>();
#>  179 :       out__.write(y);
#>  180 :       if (stan::math::logical_negation(
#>  181 :             (stan::math::primitive_value(emit_transformed_parameters__) ||
#>  182 :             stan::math::primitive_value(emit_generated_quantities__)))) {
#>  183 :         return ;
#>  184 :       }
#>  185 :       if (stan::math::logical_negation(emit_generated_quantities__)) {
#>  186 :         return ;
#>  187 :       }
#>  188 :     } catch (const std::exception& e) {
#>  189 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  190 :     }
#>  191 :   }
#>  192 :   template <typename VecVar, typename VecI,
#>  193 :             stan::require_vector_t<VecVar>* = nullptr,
#>  194 :             stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
#>  195 :   inline void
#>  196 :   unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
#>  197 :                          VecVar& vars__, std::ostream* pstream__ = nullptr) const {
#>  198 :     using local_scalar_t__ = double;
#>  199 :     stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
#>  200 :     stan::io::serializer<local_scalar_t__> out__(vars__);
#>  201 :     int current_statement__ = 0;
#>  202 :     // suppress unused var warning
#>  203 :     (void) current_statement__;
#>  204 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>  205 :     // suppress unused var warning
#>  206 :     (void) DUMMY_VAR__;
#>  207 :     try {
#>  208 :       local_scalar_t__ y = DUMMY_VAR__;
#>  209 :       current_statement__ = 1;
#>  210 :       y = in__.read<local_scalar_t__>();
#>  211 :       out__.write(y);
#>  212 :     } catch (const std::exception& e) {
#>  213 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  214 :     }
#>  215 :   }
#>  216 :   template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
#>  217 :   inline void
#>  218 :   transform_inits_impl(const stan::io::var_context& context__, VecVar&
#>  219 :                        vars__, std::ostream* pstream__ = nullptr) const {
#>  220 :     using local_scalar_t__ = double;
#>  221 :     stan::io::serializer<local_scalar_t__> out__(vars__);
#>  222 :     int current_statement__ = 0;
#>  223 :     // suppress unused var warning
#>  224 :     (void) current_statement__;
#>  225 :     local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
#>  226 :     // suppress unused var warning
#>  227 :     (void) DUMMY_VAR__;
#>  228 :     try {
#>  229 :       current_statement__ = 1;
#>  230 :       context__.validate_dims("parameter initialization", "y", "double",
#>  231 :         std::vector<size_t>{});
#>  232 :       local_scalar_t__ y = DUMMY_VAR__;
#>  233 :       current_statement__ = 1;
#>  234 :       y = context__.vals_r("y")[(1 - 1)];
#>  235 :       out__.write(y);
#>  236 :     } catch (const std::exception& e) {
#>  237 :       stan::lang::rethrow_located(e, locations_array__[current_statement__]);
#>  238 :     }
#>  239 :   }
#>  240 :   inline void
#>  241 :   get_param_names(std::vector<std::string>& names__, const bool
#>  242 :                   emit_transformed_parameters__ = true, const bool
#>  243 :                   emit_generated_quantities__ = true) const {
#>  244 :     names__ = std::vector<std::string>{"y"};
#>  245 :     if (emit_transformed_parameters__) {}
#>  246 :     if (emit_generated_quantities__) {}
#>  247 :   }
#>  248 :   inline void
#>  249 :   get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
#>  250 :            emit_transformed_parameters__ = true, const bool
#>  251 :            emit_generated_quantities__ = true) const {
#>  252 :     dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{}};
#>  253 :     if (emit_transformed_parameters__) {}
#>  254 :     if (emit_generated_quantities__) {}
#>  255 :   }
#>  256 :   inline void
#>  257 :   constrained_param_names(std::vector<std::string>& param_names__, bool
#>  258 :                           emit_transformed_parameters__ = true, bool
#>  259 :                           emit_generated_quantities__ = true) const final {
#>  260 :     param_names__.emplace_back(std::string() + "y");
#>  261 :     if (emit_transformed_parameters__) {}
#>  262 :     if (emit_generated_quantities__) {}
#>  263 :   }
#>  264 :   inline void
#>  265 :   unconstrained_param_names(std::vector<std::string>& param_names__, bool
#>  266 :                             emit_transformed_parameters__ = true, bool
#>  267 :                             emit_generated_quantities__ = true) const final {
#>  268 :     param_names__.emplace_back(std::string() + "y");
#>  269 :     if (emit_transformed_parameters__) {}
#>  270 :     if (emit_generated_quantities__) {}
#>  271 :   }
#>  272 :   inline std::string get_constrained_sizedtypes() const {
#>  273 :     return std::string("[{\"name\":\"y\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]");
#>  274 :   }
#>  275 :   inline std::string get_unconstrained_sizedtypes() const {
#>  276 :     return std::string("[{\"name\":\"y\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]");
#>  277 :   }
#>  278 :   // Begin method overload boilerplate
#>  279 :   template <typename RNG> inline void
#>  280 :   write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
#>  281 :               Eigen::Matrix<double,-1,1>& vars, const bool
#>  282 :               emit_transformed_parameters = true, const bool
#>  283 :               emit_generated_quantities = true, std::ostream*
#>  284 :               pstream = nullptr) const {
#>  285 :     const size_t num_params__ = 1;
#>  286 :     const size_t num_transformed = emit_transformed_parameters * (0U);
#>  287 :     const size_t num_gen_quantities = emit_generated_quantities * (0U);
#>  288 :     const size_t num_to_write = num_params__ + num_transformed +
#>  289 :       num_gen_quantities;
#>  290 :     std::vector<int> params_i;
#>  291 :     vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
#>  292 :              std::numeric_limits<double>::quiet_NaN());
#>  293 :     write_array_impl(base_rng, params_r, params_i, vars,
#>  294 :       emit_transformed_parameters, emit_generated_quantities, pstream);
#>  295 :   }
#>  296 :   template <typename RNG> inline void
#>  297 :   write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
#>  298 :               params_i, std::vector<double>& vars, bool
#>  299 :               emit_transformed_parameters = true, bool
#>  300 :               emit_generated_quantities = true, std::ostream*
#>  301 :               pstream = nullptr) const {
#>  302 :     const size_t num_params__ = 1;
#>  303 :     const size_t num_transformed = emit_transformed_parameters * (0U);
#>  304 :     const size_t num_gen_quantities = emit_generated_quantities * (0U);
#>  305 :     const size_t num_to_write = num_params__ + num_transformed +
#>  306 :       num_gen_quantities;
#>  307 :     vars = std::vector<double>(num_to_write,
#>  308 :              std::numeric_limits<double>::quiet_NaN());
#>  309 :     write_array_impl(base_rng, params_r, params_i, vars,
#>  310 :       emit_transformed_parameters, emit_generated_quantities, pstream);
#>  311 :   }
#>  312 :   template <bool propto__, bool jacobian__, typename T_> inline T_
#>  313 :   log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
#>  314 :     Eigen::Matrix<int,-1,1> params_i;
#>  315 :     return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
#>  316 :   }
#>  317 :   template <bool propto__, bool jacobian__, typename T_> inline T_
#>  318 :   log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
#>  319 :            std::ostream* pstream = nullptr) const {
#>  320 :     return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
#>  321 :   }
#>  322 :   inline void
#>  323 :   transform_inits(const stan::io::var_context& context,
#>  324 :                   Eigen::Matrix<double,-1,1>& params_r, std::ostream*
#>  325 :                   pstream = nullptr) const final {
#>  326 :     std::vector<double> params_r_vec(params_r.size());
#>  327 :     std::vector<int> params_i;
#>  328 :     transform_inits(context, params_i, params_r_vec, pstream);
#>  329 :     params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
#>  330 :                  params_r_vec.size());
#>  331 :   }
#>  332 :   inline void
#>  333 :   transform_inits(const stan::io::var_context& context, std::vector<int>&
#>  334 :                   params_i, std::vector<double>& vars, std::ostream*
#>  335 :                   pstream__ = nullptr) const {
#>  336 :     vars.resize(num_params_r__);
#>  337 :     transform_inits_impl(context, vars, pstream__);
#>  338 :   }
#>  339 :   inline void
#>  340 :   unconstrain_array(const std::vector<double>& params_constrained,
#>  341 :                     std::vector<double>& params_unconstrained, std::ostream*
#>  342 :                     pstream = nullptr) const {
#>  343 :     const std::vector<int> params_i;
#>  344 :     params_unconstrained = std::vector<double>(num_params_r__,
#>  345 :                              std::numeric_limits<double>::quiet_NaN());
#>  346 :     unconstrain_array_impl(params_constrained, params_i,
#>  347 :       params_unconstrained, pstream);
#>  348 :   }
#>  349 :   inline void
#>  350 :   unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
#>  351 :                     Eigen::Matrix<double,-1,1>& params_unconstrained,
#>  352 :                     std::ostream* pstream = nullptr) const {
#>  353 :     const std::vector<int> params_i;
#>  354 :     params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
#>  355 :                              std::numeric_limits<double>::quiet_NaN());
#>  356 :     unconstrain_array_impl(params_constrained, params_i,
#>  357 :       params_unconstrained, pstream);
#>  358 :   }
#>  359 : };
#>  360 : }
#>  361 : using stan_model = model1fb45d036b29__namespace::model1fb45d036b29_;
#>  362 : #ifndef USING_R
#>  363 : // Boilerplate
#>  364 : stan::model::model_base&
#>  365 : new_model(stan::io::var_context& data_context, unsigned int seed,
#>  366 :           std::ostream* msg_stream) {
#>  367 :   stan_model* m = new stan_model(data_context, seed, msg_stream);
#>  368 :   return *m;
#>  369 : }
#>  370 : stan::math::profile_map& get_stan_profile_data() {
#>  371 :   return model1fb45d036b29__namespace::profiles__;
#>  372 : }
#>  373 : #endif
#>  374 : #endif
#>  375 : 
#>  376 : RCPP_MODULE(stan_fit4model1fb45d036b29__mod) {
#>  377 :   class_<rstan::stan_fit<stan_model, boost::random::mixmax> >(
#>  378 :       "stan_fit4model1fb45d036b29_")
#>  379 : 
#>  380 :       .constructor<SEXP, SEXP, SEXP>()
#>  381 : 
#>  382 :       .method(
#>  383 :           "call_sampler",
#>  384 :           &rstan::stan_fit<stan_model, boost::random::mixmax>::call_sampler)
#>  385 :       .method(
#>  386 :           "param_names",
#>  387 :           &rstan::stan_fit<stan_model, boost::random::mixmax>::param_names)
#>  388 :       .method("param_names_oi",
#>  389 :               &rstan::stan_fit<stan_model,
#>  390 :                                boost::random::mixmax>::param_names_oi)
#>  391 :       .method("param_fnames_oi",
#>  392 :               &rstan::stan_fit<stan_model,
#>  393 :                                boost::random::mixmax>::param_fnames_oi)
#>  394 :       .method(
#>  395 :           "param_dims",
#>  396 :           &rstan::stan_fit<stan_model, boost::random::mixmax>::param_dims)
#>  397 :       .method("param_dims_oi",
#>  398 :               &rstan::stan_fit<stan_model,
#>  399 :                                boost::random::mixmax>::param_dims_oi)
#>  400 :       .method("update_param_oi",
#>  401 :               &rstan::stan_fit<stan_model,
#>  402 :                                boost::random::mixmax>::update_param_oi)
#>  403 :       .method("param_oi_tidx",
#>  404 :               &rstan::stan_fit<stan_model,
#>  405 :                                boost::random::mixmax>::param_oi_tidx)
#>  406 :       .method("grad_log_prob",
#>  407 :               &rstan::stan_fit<stan_model,
#>  408 :                                boost::random::mixmax>::grad_log_prob)
#>  409 :       .method("log_prob",
#>  410 :               &rstan::stan_fit<stan_model, boost::random::mixmax>::log_prob)
#>  411 :       .method("unconstrain_pars",
#>  412 :               &rstan::stan_fit<stan_model,
#>  413 :                                boost::random::mixmax>::unconstrain_pars)
#>  414 :       .method("constrain_pars",
#>  415 :               &rstan::stan_fit<stan_model,
#>  416 :                                boost::random::mixmax>::constrain_pars)
#>  417 :       .method(
#>  418 :           "num_pars_unconstrained",
#>  419 :           &rstan::stan_fit<stan_model,
#>  420 :                            boost::random::mixmax>::num_pars_unconstrained)
#>  421 :       .method(
#>  422 :           "unconstrained_param_names",
#>  423 :           &rstan::stan_fit<
#>  424 :               stan_model, boost::random::mixmax>::unconstrained_param_names)
#>  425 :       .method(
#>  426 :           "constrained_param_names",
#>  427 :           &rstan::stan_fit<stan_model,
#>  428 :                            boost::random::mixmax>::constrained_param_names)
#>  429 :       .method("standalone_gqs",
#>  430 :               &rstan::stan_fit<stan_model,
#>  431 :                                boost::random::mixmax>::standalone_gqs);
#>  432 : }
#>  433 : 
#>  434 : 
#>  435 : // declarations
#>  436 : extern "C" {
#>  437 : SEXP file1fb4623e2561( ) ;
#>  438 : }
#>  439 : 
#>  440 : // definition
#>  441 : SEXP file1fb4623e2561() {
#>  442 :  return Rcpp::wrap("anon_model");
#>  443 : }
fit <- sampling(mod, data = list(y_mean = 0))
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 1e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.01 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.003 seconds (Warm-up)
#> Chain 1:                0.002 seconds (Sampling)
#> Chain 1:                0.005 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.003 seconds (Warm-up)
#> Chain 2:                0.003 seconds (Sampling)
#> Chain 2:                0.006 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.003 seconds (Warm-up)
#> Chain 3:                0.003 seconds (Sampling)
#> Chain 3:                0.006 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 0 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.003 seconds (Warm-up)
#> Chain 4:                0.003 seconds (Sampling)
#> Chain 4:                0.006 seconds (Total)
#> Chain 4: 
fit2 <- sampling(mod, data = list(y_mean = 5))
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 0 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.003 seconds (Warm-up)
#> Chain 1:                0.003 seconds (Sampling)
#> Chain 1:                0.006 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 0 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.003 seconds (Warm-up)
#> Chain 2:                0.003 seconds (Sampling)
#> Chain 2:                0.006 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 0 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.003 seconds (Warm-up)
#> Chain 3:                0.003 seconds (Sampling)
#> Chain 3:                0.006 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'anon_model' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 0 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.003 seconds (Warm-up)
#> Chain 4:                0.003 seconds (Sampling)
#> Chain 4:                0.006 seconds (Total)
#> Chain 4: 
# }
```
