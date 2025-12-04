# Run Stan's variational algorithm for approximate posterior sampling

Approximately draw from a posterior distribution using variational
inference.

This is still considered an experimental feature. We recommend calling
[`stan`](https://mc-stan.org/rstan/dev/reference/stan.md) or
[`sampling`](https://mc-stan.org/rstan/dev/reference/stanmodel-method-sampling.md)
for final inferences and only using `vb` to get a rough idea of the
parameter distributions.

## Usage

    <!-- %% vb(object, \dots)   -->
      # S4 method for class 'stanmodel'
    vb(object, data = list(), pars = NA, include = TRUE,
        seed = sample.int(.Machine$integer.max, 1),
        init = 'random', check_data = TRUE,
        sample_file = tempfile(fileext = '.csv'),
        algorithm = c("meanfield", "fullrank"),
        importance_resampling = FALSE, keep_every = 1,
        ...)

## Methods

- vb:

  `signature(object = "stanmodel")` Call Stan's variational Bayes
  methods for the model defined by S4 class `stanmodel` given the data,
  initial values, etc.

## Arguments

- object:

  An object of class
  [`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md).

- data:

  A named `list` or `environment` providing the data for the model or a
  character vector for all the names of objects used as data. See the
  **Passing data to Stan** section in
  [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md).

- pars:

  If not `NA`, then a character vector naming parameters, which are
  included in the output if `include = TRUE` and excluded if
  `include = FALSE`. By default, all parameters are included.

- include:

  Logical scalar defaulting to `TRUE` indicating whether to include or
  exclude the parameters given by the `pars` argument. If `FALSE`, only
  entire multidimensional parameters can be excluded, rather than
  particular elements of them.

- seed:

  The seed for random number generation. The default is generated from 1
  to the maximum integer supported by R on the machine. Even if multiple
  chains are used, only one seed is needed, with other chains having
  seeds derived from that of the first chain to avoid dependent samples.
  When a seed is specified by a number, `as.integer` will be applied to
  it. If `as.integer` produces `NA`, the seed is generated randomly. The
  seed can also be specified as a character string of digits, such as
  `"12345"`, which is converted to integer.

- init:

  Initial values specification. See the detailed documentation for the
  init argument in
  [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md).

- check_data:

  Logical, defaulting to `TRUE`. If `TRUE` the data will be
  preprocessed; otherwise not. See the **Passing data to Stan** section
  in [`stan`](https://mc-stan.org/rstan/dev/reference/stan.md).

- sample_file:

  A character string of file name for specifying where to write samples
  for *all* parameters and other saved quantities. This defaults to a
  temporary file.

- algorithm:

  Either `"meanfield"` (the default) or `"fullrank"`, indicating which
  variational inference algorithm is used. The `"meanfield"` option uses
  a fully factorized Gaussian for the approximation whereas the
  `fullrank` option uses a Gaussian with a full-rank covariance matrix
  for the approximation. Details and additional references are available
  in the Stan manual.

- importance_resampling:

  Logical scalar (defaulting to `FALSE`) indicating whether to do
  importance resampling to adjust the draws at the optimum to be more
  like draws from the posterior distribution

- keep_every:

  Integer scalar (defaulting to 1) indicating the interval by which to
  thin the draws when `imporance_resampling = TRUE`

- ...:

  Other optional parameters:

  - `iter` (positive `integer`), the maximum number of iterations,
    defaulting to 10000.

  - `grad_samples` (positive `integer`), the number of samples for Monte
    Carlo estimate of gradients, defaulting to 1.

  - `elbo_samples` (positive `integer`), the number of samples for Monte
    Carlo estimate of ELBO (objective function), defaulting to 100.
    (ELBO stands for "the evidence lower bound".)

  - `eta` (`double`), positive stepsize weighting parameter for
    variational inference but is ignored if adaptation is engaged, which
    is the case by default.

  - `adapt_engaged` (`logical`), a flag indicating whether to
    automatically adapt the stepsize, defaulting to `TRUE`.

  - `tol_rel_obj` (positive `double`), the convergence tolerance on the
    relative norm of the objective, defaulting to 0.01.

  - `eval_elbo` (positive `integer`), evaluate ELBO every Nth iteration,
    defaulting to 100.

  - `output_samples` (positive `integer`), number of posterior samples
    to draw and save, defaults to 1000.

  - `adapt_iter` (positive `integer`), the maximum number of iterations
    to adapt the stepsize, defaulting to 50. Ignored if
    `adapt_engaged = FALSE`.

  Refer to the manuals for both CmdStan and Stan for more details.

## Value

An object of
[`stanfit-class`](https://mc-stan.org/rstan/dev/reference/stanfit-class.md).

## See also

[`stanmodel`](https://mc-stan.org/rstan/dev/reference/stanmodel-class.md)

The manuals of CmdStan and Stan.

## References

The Stan Development Team *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org>.

The Stan Development Team *CmdStan Interface User's Guide*.
<https://mc-stan.org>.

## Examples

``` r
# \dontrun{
m <- stan_model(model_code = 'parameters {real y;} model {y ~ normal(0,1);}')
#> recompiling to avoid crashing R session
f <- vb(m)
#> Chain 1: ------------------------------------------------------------
#> Chain 1: EXPERIMENTAL ALGORITHM:
#> Chain 1:   This procedure has not been thoroughly tested and may be unstable
#> Chain 1:   or buggy. The interface is subject to change.
#> Chain 1: ------------------------------------------------------------
#> Chain 1: 
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2e-06 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.02 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Begin eta adaptation.
#> Chain 1: Iteration:   1 / 250 [  0%]  (Adaptation)
#> Chain 1: Iteration:  50 / 250 [ 20%]  (Adaptation)
#> Chain 1: Iteration: 100 / 250 [ 40%]  (Adaptation)
#> Chain 1: Iteration: 150 / 250 [ 60%]  (Adaptation)
#> Chain 1: Iteration: 200 / 250 [ 80%]  (Adaptation)
#> Chain 1: Success! Found best value [eta = 1] earlier than expected.
#> Chain 1: 
#> Chain 1: Begin stochastic gradient ascent.
#> Chain 1:   iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes 
#> Chain 1:    100           -0.033             1.000            1.000
#> Chain 1:    200            0.078             1.213            1.425
#> Chain 1:    300            0.104             0.892            1.000
#> Chain 1:    400           -0.092             1.200            1.425
#> Chain 1:    500            0.052             1.517            1.425
#> Chain 1:    600            0.057             1.279            1.425
#> Chain 1:    700           -0.072             1.352            1.425
#> Chain 1:    800           -0.011             1.896            1.791
#> Chain 1:    900           -0.038             1.765            1.425
#> Chain 1:   1000            0.039             1.786            1.791
#> Chain 1:   1100            0.077             1.736            1.791   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1200           -0.238             1.726            1.791   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1300            0.021             2.921            1.976   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1400            0.088             2.785            1.791   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1500           -0.162             2.661            1.547   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1600           -0.060             2.819            1.671   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1700           -0.114             2.687            1.547   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1800            0.100             2.331            1.547   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   1900            0.089             2.272            1.547   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2000            0.102             2.087            1.323   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2100           -0.046             2.358            1.547   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2200           -0.068             2.257            1.547   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2300            0.008             1.960            1.547   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2400           -0.165             1.990            1.547   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2500           -0.007             4.246            1.671   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2600           -0.014             4.132            1.050   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2700            0.083             4.202            1.170   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2800           -0.020             4.494            1.170   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   2900            0.037             4.636            1.544   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3000            0.013             4.801            1.778   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3100            0.064             4.559            1.544   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3200            0.065             4.528            1.544   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3300            0.022             3.800            1.544   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3400            0.094             3.772            1.544   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3500            0.046             1.464            1.170   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3600            0.014             1.643            1.544   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3700           -0.041             1.660            1.544   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3800           -0.006             1.768            1.544   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   3900           -0.030             1.695            1.338   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4000           -0.084             1.581            1.034   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4100           -0.011             2.192            1.338   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4200           -0.000             7.171            1.959   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4300           -0.044             7.075            1.338   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4400            0.135             7.131            1.338   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4500           -0.059             7.357            2.325   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4600           -0.139             7.182            1.338   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4700            0.014             8.154            3.292   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4800            0.051             7.613            1.327   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   4900           -0.020             7.891            3.292   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5000           -0.034             7.870            3.292   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5100            0.075             7.325            1.452   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5200            0.099             2.369            1.327   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5300           -0.138             2.441            1.452   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5400           -0.102             2.344            1.452   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5500           -0.234             2.071            0.728   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5600            0.115             2.316            1.452   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5700           -0.016             2.047            1.452   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5800           -0.141             2.063            1.452   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   5900           -0.003             6.437            1.452   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6000           -0.019             6.479            1.452   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6100           -0.068             6.406            0.889   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6200           -0.007             7.231            1.714   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6300           -0.061             7.148            0.889   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6400           -0.041             7.164            0.889   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6500           -0.043             7.114            0.889   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6600           -0.003             8.161            0.889   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6700            0.084             7.428            0.889   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6800            0.001            17.772            1.035   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   6900            0.054            13.137            0.985   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7000           -0.044            13.274            1.035   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7100            0.030            13.450            2.219   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7200           -0.035            12.786            1.845   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7300           -0.042            12.714            1.845   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7400           -0.110            12.724            1.845   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7500            0.069            12.977            2.219   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7600           -0.101            11.796            1.845   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7700            0.119            11.877            1.845   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7800            0.001            22.918            1.845   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   7900           -0.005            22.929            1.845   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8000            0.024            22.830            1.845   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8100           -0.137            22.699            1.689   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8200            0.046            22.912            1.689   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8300            0.044            22.901            1.689   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8400           -0.036            23.059            1.845   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8500           -0.005            23.373            1.845   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8600            0.117            23.309            1.845   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8700            0.038            23.333            2.087   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8800           -0.132             1.988            1.289   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   8900           -0.119             1.889            1.289   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9000            0.052             2.094            2.087   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9100            0.050             1.980            2.087   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9200           -0.031             1.845            2.087   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9300            0.025             2.063            2.194   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9400           -0.073             1.978            2.087   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9500           -0.060             1.426            1.344   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9600            0.085             1.492            1.706   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9700            0.097             1.296            1.344   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9800            0.054             1.247            1.344   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   9900            0.030             1.316            1.344   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1:   10000           -0.009             1.413            1.344   MAY BE DIVERGING... INSPECT ELBO
#> Chain 1: Informational Message: The maximum number of iterations is reached! The algorithm may not have converged.
#> Chain 1: This variational approximation is not guaranteed to be meaningful.
#> Chain 1: 
#> Chain 1: Drawing a sample of size 1000 from the approximate posterior... 
#> Chain 1: COMPLETED.
# }
```
