files <- Sys.glob("../lib/libStanHeaders.a")
stopifnot(length(files) > 0)
dest <- file.path(R_PACKAGE_DIR, paste0('lib', R_ARCH))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
file.copy(files, dest, overwrite = TRUE)
if (file.exists("symbols.rds"))
  file.copy("symbols.rds", dest, overwrite = TRUE)

