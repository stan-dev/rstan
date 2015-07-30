# This file is part of RStan
# Copyright (C) 2012, 2013 Hadley Wickham
# Copyright (C) 2015 Jiqiang Guo and Benjamin Goodrich
#
# RStan is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# RStan is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

find_rtools <- function(debug = FALSE) {
  # Non-windows users don't need rtools
  if (.Platform$OS.type != "windows") return(TRUE)
  
  
  # First try the path
  from_path <- scan_path_for_rtools(debug)
  if (is_compatible(from_path)) return(TRUE)

  if (!is.null(from_path)) {
    # Installed
    if (is.null(from_path$version)) {
      # but not from rtools
      if (debug) "gcc and ls on path, assuming set up is correct\n"
      return(TRUE)
    } else {
      # Installed, but not compatible
      warning("WARNING: Rtools ", from_path$version, " found on the path",
              " at ", from_path$path, " is not compatible with R ", getRversion(), ".\n\n",
              "Please download and install ", rtools_needed(), " from ", rtools_url,
              ", remove the incompatible version from your PATH, then run find_rtools().")
      return(invisible(FALSE))
    }
  }
  
  # Not on path, so try registry
  registry_candidates <- scan_registry_for_rtools(debug)
  
  if (length(registry_candidates) == 0) {
    # Not on path or in registry, so not installled
    warning("WARNING: Rtools is required to build Stan programs, but is not ",
            "currently installed.\n\n",
            "Please download and install ", rtools_needed(), " from ", rtools_url, ".")
    return(invisible(FALSE))
  }
  
  from_registry <- Find(is_compatible, registry_candidates, right = TRUE)
  if (is.null(from_registry)) {
    # In registry, but not compatible.
    versions <- vapply(registry_candidates, function(x) x$version, character(1))
    warning("WARNING: Rtools is required to build Stan programs, but no version ",
            "of Rtools compatible with R ", getRversion(), " was found. ",
            "(Only the following incompatible version(s) of Rtools were found:",
            paste(versions, collapse = ","), ")\n\n",
            "Please download and install ", rtools_needed(), " from ", rtools_url, ".")
    return(invisible(FALSE))
  }
  
  installed_ver <- installed_version(from_registry$path, debug = debug)
  if (is.null(installed_ver)) {
    # Previously installed version now deleted
    warning("WARNING: Rtools is required to build Stan programs, but the ",
            "version of Rtools previously installed in ", from_registry$path,
            " has been deleted.\n\n",
            "Please download and install ", rtools_needed(), " from ", rtools_url, ".")
    return(invisible(FALSE))
  }
  
  if (installed_ver != from_registry$version) {
    # Installed version doesn't match registry version
    warning("WARNING: Rtools is required to build Stan programs, but no version ",
            "of Rtools compatible with R ", getRversion(), " was found. ",
            "Rtools ", from_registry$version, " was previously installed in ",
            from_registry$path, " but now that directory contains Rtools ",
            installed_ver, ".\n\n",
            "Please download and install ", rtools_needed(), " from ", rtools_url, ".")
    return(invisible(FALSE))
  }
  
  TRUE
}

scan_path_for_rtools <- function(debug = FALSE) {
  if (debug) cat("Scanning path...\n")
  
  # First look for ls and gcc
  ls_path <- Sys.which("ls")
  if (ls_path == "") return(NULL)
  if (debug) cat("ls :", ls_path, "\n")
  
  gcc_path <- Sys.which("gcc")
  if (gcc_path == "") return(NULL)
  if (debug) cat("gcc:", gcc_path, "\n")
  
  # We have a candidate installPath
  install_path <- dirname(dirname(ls_path))
  install_path2 <- dirname(dirname(dirname(gcc_path)))
  if (tolower(install_path2) != tolower(install_path)) return(NULL)
  
  version <- installed_version(install_path, debug = debug)
  if (debug) cat("Version:", version, "\n")
  
  rtools(install_path, version)
}

scan_registry_for_rtools <- function(debug = FALSE) {
  if (debug) cat("Scanning registry...\n")
  
  keys <- NULL
  try(keys <- utils::readRegistry("SOFTWARE\\R-core\\Rtools",
                                  hive = "HCU", view = "32-bit", maxdepth = 2), silent = TRUE)
  if (is.null(keys))
    try(keys <- utils::readRegistry("SOFTWARE\\R-core\\Rtools",
                                    hive = "HLM", view = "32-bit", maxdepth = 2), silent = TRUE)
  if (is.null(keys)) return(NULL)
  
  rts <- vector("list", length(keys))
  
  for(i in seq_along(keys)) {
    version <- names(keys)[[i]]
    key <- keys[[version]]
    if (!is.list(key) || is.null(key$InstallPath)) next;
    install_path <- normalizePath(key$InstallPath,
                                  mustWork = FALSE, winslash = "/")
    
    if (debug) cat("Found", install_path, "for", version, "\n")
    rts[[i]] <- rtools(install_path, version)
  }
  
  Filter(Negate(is.null), rts)
}

installed_version <- function(path, debug) {
  if (!file.exists(file.path(path, "Rtools.txt"))) return(NULL)
  
  # Find the version path
  version_path <- file.path(path, "VERSION.txt")
  if (debug) {
    cat("VERSION.txt\n")
    cat(readLines(version_path), "\n")
  }
  if (!file.exists(version_path)) return(NULL)
  
  # Rtools is in the path -- now crack the VERSION file
  contents <- NULL
  try(contents <- readLines(version_path), silent = TRUE)
  if (is.null(contents)) return(NULL)
  
  # Extract the version
  contents <- gsub("^\\s+|\\s+$", "", contents)
  version_re <- "Rtools version (\\d\\.\\d+)\\.[0-9.]+$"
  
  if (!grepl(version_re, contents)) return(NULL)
  
  m <- regexec(version_re, contents)
  regmatches(contents, m)[[1]][2]
}

is_compatible <- function(rtools) {
  if (is.null(rtools)) return(FALSE)
  if (is.null(rtools$version)) return(FALSE)
  
  stopifnot(is.rtools(rtools))
  info <- version_info[[rtools$version]]
  if (is.null(info)) return(FALSE)
  
  r_version <- getRversion()
  r_version >= info$version_min && r_version <= info$version_max
}

rtools <- function(path, version) {
  structure(list(version = version, path = path), class = "rtools")
}
is.rtools <- function(x) inherits(x, "rtools")

# Rtools metadata --------------------------------------------------------------
rtools_url <- "http://cran.r-project.org/bin/windows/Rtools/"
version_info <- list(
  "2.11" = list(
    version_min = "2.10.0",
    version_max = "2.11.1",
    path = c("bin", "perl/bin", "MinGW/bin")
  ),
  "2.12" = list(
    version_min = "2.12.0",
    version_max = "2.12.2",
    path = c("bin", "perl/bin", "MinGW/bin", "MinGW64/bin")
  ),
  "2.13" = list(
    version_min = "2.13.0",
    version_max = "2.13.2",
    path = c("bin", "MinGW/bin", "MinGW64/bin")
  ),
  "2.14" = list(
    version_min = "2.13.0",
    version_max = "2.14.2",
    path = c("bin", "MinGW/bin", "MinGW64/bin")
  ),
  "2.15" = list(
    version_min = "2.14.2",
    version_max = "2.15.1",
    path = c("bin", "gcc-4.6.3/bin")
  ),
  "2.16" = list(
    version_min = "2.15.2",
    version_max = "3.0.0",
    path = c("bin", "gcc-4.6.3/bin")
  ),
  "3.0" = list(
    version_min = "2.15.2",
    version_max = "3.0.99",
    path = c("bin", "gcc-4.6.3/bin")
  ),
  "3.1" = list(
    version_min = "3.0.0",
    version_max = "3.1.99",
    path = c("bin", "gcc-4.6.3/bin")
  ),
  "3.2" = list(
    version_min = "3.1.0",
    version_max = "3.2.99",
    path = c("bin", "gcc-4.6.3/bin")
  ),
  "3.3" = list(
    version_min = "3.2.0",
    version_max = "3.3.99",
    path = c("bin", "gcc-4.6.3/bin")
  )
)

rtools_needed <- function() {
  r_version <- getRversion()
  
  for(i in rev(seq_along(version_info))) {
    version <- names(version_info)[i]
    info <- version_info[[i]]
    ok <- r_version >= info$version_min && r_version <= info$version_max
    if (ok) return(paste("Rtools", version))
  }
  "the appropriate version of Rtools"
}
