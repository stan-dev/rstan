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



# defaults ----------------------------------------------------------------
.rstanvis_defaults <- new.env(parent = emptyenv())
.rstanvis_defaults$theme <-
  theme_classic(base_size = 11, base_family = "Arial") +
  theme(axis.line = element_line(color = "#222222"),
        axis.line.y = element_line(size = .5),
        axis.line.x = element_line(size = 1),
        axis.title = element_text(face = "bold", size = 13),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", face = "bold"),
        legend.title = element_text(size = 11),
        plot.title = element_text(size = 18))

for (j in seq_along(.rstanvis_defaults$theme)) {
  if ("element_text" %in% class(.rstanvis_defaults$theme[[j]])) {
    .rstanvis_defaults$theme[[j]][["debug"]] <- FALSE
  }
}

.rstanvis_defaults$hist_theme <-
  .rstanvis_defaults$theme +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank())

.rstanvis_defaults$multiparam_theme <-
  .rstanvis_defaults$theme +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(face = "bold", size = 13, debug = FALSE),
        legend.position = "none",
        panel.grid.major =  element_line(size = 0.1, color = "darkgray"))

.rstanvis_defaults$pt_color <- "#B2001D"
.rstanvis_defaults$alpha <- 0.5
.rstanvis_defaults$shape <- 19
.rstanvis_defaults$fill <-  "#B2001D"
.rstanvis_defaults$color <- "black" #"#590815"
.rstanvis_defaults$size <- 0.5
.rstanvis_defaults$pt_size <- 3

.rstanvis_defaults$chain_colors <- rgb(matrix(c(230, 97, 1,
                                                153, 142, 195,
                                                84, 39, 136,
                                                241, 163, 64,
                                                216, 218, 235,
                                                254, 224, 182),
                                              byrow = TRUE, ncol = 3),
                                       names = paste(1:6), maxColorValue = 255)
.rstanvis_defaults$grays <- rgb(matrix(c(247, 247, 247, 204, 204, 204,
                                         150, 150, 150, 82, 82, 82),
                                       byrow = TRUE, ncol = 3),
                                alpha = 100, names = paste(1:4), maxColorValue = 255)
