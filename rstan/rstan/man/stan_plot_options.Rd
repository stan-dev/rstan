\name{rstan_gg_options}
\alias{rstan_gg_options}
\alias{rstan_ggtheme_options}
\title{Set default appearance options}
\usage{
  rstan_gg_options(...)
  
  rstan_ggtheme_options(...)
}
\arguments{
  \item{...}{For \code{rstan_ggtheme_options}, see
    \code{\link[ggplot2]{theme}} for the theme elements that can be specified
    in \code{...}. For \code{rstan_gg_options}, \code{...} can be \code{fill},
    \code{color}, \code{chain_colors}, or \code{size}. See Examples.}
}
\description{
  Set default appearance options
}
\seealso{
 \code{\link[=stan_plot]{Plots}}, \code{\link[=stan_diag]{Diagnostic Plots}}
}
\examples{

rstan_ggtheme_options(panel.background = element_rect(fill = "gray"),
                      legend.position = "top")
rstan_gg_options(fill = "skyblue", color = "skyblue4")
}

