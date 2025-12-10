# Set default appearance options

Set default appearance options

## Usage

``` r
rstan_gg_options(...)
  
  rstan_ggtheme_options(...)
```

## Arguments

- ...:

  For `rstan_ggtheme_options`, see
  [`theme`](https://ggplot2.tidyverse.org/reference/theme.html) for the
  theme elements that can be specified in `...`. For `rstan_gg_options`,
  `...` can be `fill`, `color`, `chain_colors`, `size`, `pt_color`, or
  `pt_size`. See Examples.

## See also

[`List of RStan plotting functions`](https://mc-stan.org/rstan/reference/plotting-functions.md)

## Examples

``` r
rstan_ggtheme_options(panel.background = ggplot2::element_rect(fill = "gray"),
                      legend.position = "top")
rstan_gg_options(fill = "skyblue", color = "skyblue4", pt_color = "red")
```
