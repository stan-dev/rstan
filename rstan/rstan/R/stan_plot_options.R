rstan_gg_options <- function(...) {
  ops <- c("size", "fill", "color", "chain_colors")
  dots <- list(...)
  if (length(dots) == 0) {
    oplist <- list()
    for (op in ops) oplist[[op]] <- .rstanvis_defaults[[op]]
    print(oplist)
    return(invisible(NULL))
  }
  nms <- names(dots)
  for (x in nms) {
    .rstanvis_defaults[[x]] <- dots[[x]]
  }
}

rstan_ggtheme_options <- function(...) {
  if (length(list(...)) == 0) {
    print(.rstanvis_defaults$theme)
    return(invisible(NULL))
  }
  els <- theme(...)
  nms <- names(els)
  for (x in nms) {
    .rstanvis_defaults$theme[[x]] <- els[[x]]
    .rstanvis_defaults$hist_theme[[x]] <- els[[x]]
    .rstanvis_defaults$multiparam_theme[[x]] <- els[[x]]
  }
}



# internal ----------------------------------------------------------------

.rstanvis_defaults <- new.env(parent = emptyenv())

rstanvis_theme <- function() .rstanvis_defaults$theme
rstanvis_hist_theme <- function() .rstanvis_defaults$hist_theme
rstanvis_multiparam_theme <- function() .rstanvis_defaults$multiparam_theme
rstanvis_aes_ops <- function(nm = NULL) {
  if (is.null(nm)) return(.rstanvis_defaults$aes_ops)
  else return(.rstanvis_defaults$aes_ops[[nm]])
}

# called in .onLoad()
set_rstan_ggplot_defaults <- function() {
  .rstanvis_defaults$theme <- ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.line = ggplot2::element_line(color = "#222222"),
          axis.line.y = ggplot2::element_line(size = .5),
          axis.line.x = ggplot2::element_line(size = 1),
          axis.title = ggplot2::element_text(face = "bold", size = 13),
          strip.background = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(color = "black", face = "bold"),
          legend.title = ggplot2::element_text(size = 11),
          plot.title = ggplot2::element_text(size = 18))
  
  for (j in seq_along(.rstanvis_defaults$theme)) {
    if ("element_text" %in% class(.rstanvis_defaults$theme[[j]])) {
      .rstanvis_defaults$theme[[j]][["debug"]] <- FALSE
    }
  }
  
  .rstanvis_defaults$hist_theme <-
    .rstanvis_defaults$theme +
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.line.y = ggplot2::element_blank())
  
  .rstanvis_defaults$multiparam_theme <-
    .rstanvis_defaults$theme +
    theme(axis.title = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(face = "bold", size = 13, debug = FALSE),
          legend.position = "none",
          panel.grid.major =  ggplot2::element_line(size = 0.1, color = "darkgray"))
  
  .rstanvis_defaults$aes_ops <- list()
  .rstanvis_defaults$aes_ops$pt_color <- "#B2001D"
  .rstanvis_defaults$aes_ops$alpha <- 0.5
  .rstanvis_defaults$aes_ops$shape <- 19
  .rstanvis_defaults$aes_ops$fill <-  "#B2001D"
  .rstanvis_defaults$aes_ops$color <- "black" #"#590815"
  .rstanvis_defaults$aes_ops$size <- 0.5
  .rstanvis_defaults$aes_ops$pt_size <- 3
  .rstanvis_defaults$aes_ops$chain_colors <- 
    rgb(
      matrix(
        c(230, 97, 1,
          153, 142, 195,
          84, 39, 136,
          241, 163, 64,
          216, 218, 235,
          254, 224, 182),
        byrow = TRUE, 
        ncol = 3
      ),
      names = paste(1:6), 
      maxColorValue = 255
    )
  
  .rstanvis_defaults$aes_ops$grays <-
    rgb(
      matrix(
        c(247, 247, 247, 204, 204, 204,
          150, 150, 150, 82, 82, 82),
        byrow = TRUE,
        ncol = 3
      ),
      alpha = 100,
      names = paste(1:4),
      maxColorValue = 255
    )
}
