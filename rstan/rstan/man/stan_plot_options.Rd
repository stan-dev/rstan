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
    \code{color}, \code{chain_colors}, \code{size}, \code{pt_color}, or
    \code{pt_size}. See Examples.}
}
\description{
  Set default appearance options
}
\seealso{
\code{\link[=rstan-plotting-functions]{List of RStan plotting functions}}
}
\examples{
rstan_ggtheme_options(panel.background = ggplot2::element_rect(fill = "gray"),
                      legend.position = "top")
rstan_gg_options(fill = "skyblue", color = "skyblue4", pt_color = "red")
}

