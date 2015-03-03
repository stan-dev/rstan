stan_demo <-
function(model = character(0), 
         method = c("sampling", "optimizing"), ...) {
  if(is.numeric(model)) {
    MODEL_NUM <- as.integer(model)
    model <- character(0)
  }
  else MODEL_NUM <- -1L
  MODELS_HOME <- system.file('include', 'example-models-master', package = 'rstan')
  if(MODELS_HOME == "") {
    MODELS_HOME <- file.path(tempdir(), "example-models-master")
    if(!file.exists(MODELS_HOME)) {
      writable <- file.access(system.file('include', package = 'rstan'), 
                              mode = 2) == 0
      choices <- c(if(writable) "Download to include directory of the rstan package",
                   "Download to the temporary directory", "Do not download")
      if(interactive()) choice <- menu(choices, title = "Do you want to download the example models?")
      else choice <- if(writable) 2L else 1L
      if(choice == 0 | choice == length(choices)) stop("cannot proceed without example models")
      if(choices[choice] == "Download to include directory of the rstan package") {
        FILE <- file.path(system.file('include', package = 'rstan'),
                          "example-models-master.zip")
      }
      else FILE <- file.path(tempdir(), "example-models-master.zip")
      request <- httr::GET("https://github.com/stan-dev/example-models/archive/master.zip")
      writeBin(httr::content(request, type = "raw"), FILE)
      unzip(FILE, exdir = dirname(FILE))
      MODELS_HOME <- file.path(dirname(FILE), "example-models-master")
    }
  }
  WD <- getwd()
  on.exit(setwd(WD))
  setwd(MODELS_HOME)
  MODELS <- dir(MODELS_HOME, pattern = paste0(model, ".stan", "$"), 
                recursive = TRUE, full.names = FALSE)
  if(length(MODELS) == 0) {
    stop("'model' not found; leave 'model' unspecified to see all choices")
  }
  else if(length(MODELS) > 1) {
    if(MODEL_NUM %in%  1:length(MODELS)) {
      MODELS <- MODELS[MODEL_NUM]
    }
    else if(MODEL_NUM == 0) MODELS <- ""
    else MODELS <- select.list(MODELS)
    if(!nzchar(MODELS)) {
      return(dir(MODELS_HOME, pattern = paste0(model, ".stan", "$"), 
                 recursive = TRUE, full.names = FALSE))
    }
    model <- sub(".stan$", "", basename(MODELS))
  }
  MODEL_HOME <- dirname(MODELS)
  STAN_ENV <- new.env()
  if(file.exists(fp <- file.path(MODEL_HOME, paste0(model, ".data.R")))) {
    source(fp, local = STAN_ENV, verbose = FALSE, echo = TRUE)
  }
  method <- match.arg(method)
  dots <- list(...)
  if(is.null(dots$object)) dots$object <- stan_model(MODELS, model_name = model)
  dots$data <- STAN_ENV
  if(method == "sampling") return(do.call(sampling, args = dots))
  else return(do.call(optimizing, args = dots))
}
