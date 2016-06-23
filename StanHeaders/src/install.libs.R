if (.Platform$OS.type == "unix") {
  files <- Sys.glob("../lib/libStanHeaders.a")  
  dest <- file.path(R_PACKAGE_DIR, paste0('lib', R_ARCH))
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(files, dest, overwrite = TRUE)
  if (file.exists("symbols.rds"))
    file.copy("symbols.rds", dest, overwrite = TRUE)
} else {
  files <- Sys.glob(paste0("*", SHLIB_EXT))
  dest <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(files, dest, overwrite = TRUE)
  if(file.exists("symbols.rds"))
    file.copy("symbols.rds", dest, overwrite = TRUE)  
}

