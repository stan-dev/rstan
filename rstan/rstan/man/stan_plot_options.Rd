\name{rstanvis_options}
\alias{rstanvis_options}
\alias{rstanvis_theme_options}
\title{Set default appearance options for ggplot graphics}
\usage{
  rstanvis_options(...)
  
  rstanvis_theme_options(...)
}
\arguments{
  \item{...}{For \code{rstanvis_theme_options}, see
    \code{\link[ggplot2]{theme}} for the theme elements that can be specified
    in \code{...}. For \code{rstanvis_options}, \code{...} can be \code{fill},
    \code{color}, \code{chain_colors}, or \code{size}. See Examples.}
}
\description{
  Set Default Appearance Options
}
\examples{
  rstanvis_theme_options(panel.background = element_rect(fill = "gray"),
                         legend.position = "top")
  rstanvis_options(fill = "skyblue", color = "skyblue4")
}

