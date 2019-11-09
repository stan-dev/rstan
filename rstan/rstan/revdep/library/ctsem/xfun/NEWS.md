# CHANGES IN xfun VERSION 0.6

## NEW FEATURES

- Added the `...` argument to `same_path()` to pass additional arguments to `normalize_path()`.

## BUG FIXES

- The `warn` argument in `prose_index()` failed to suppress warnings.

# CHANGES IN xfun VERSION 0.5

## NEW FEATURES

- Added functions `upload_ftp()` and `upload_win_builder()` to upload files to FTP servers.

- Added a function `stringsAsStrings()` (see its help page for details).

- Added an argument `warn` to `prose_index()` to suppress the warning when code fences are not balanced.

## BUG FIXES

- Fixed the bug that `prose_index()` recognizes double backticks as code fences (thanks, @shrektan, #14 #15).

# CHANGES IN xfun VERSION 0.4

## NEW FEATURES

- Added functions `embed_file()`, `embed_dir()`, and `embed_files()` to embed files in an HTML output file (e.g., from R Markdown), so that the files can be directly downloaded from the web browser. One use case is to call one of these functions in an R code chunk of an Rmd document to embed the Rmd source document or data files in the HTML output, so readers can download them.

- Added a new argument `message` to `pkg_attach()`, so you can suppress package startup messages via `xfun::pkg_attach(..., message = FALSE)` or set the global option `options(xfun.pkg_attach.message = FALSE)` (thanks, @wch, yihui/knitr#1583).

## MINOR CHANGES

- The argument `rw_error` was moved from `gsub_dir()` to `gsub_file()` (`gsub_dir(rw_error = ...)` will still work).

- `is_ascii()` now returns `NA` for `NA_character_` (thanks, @shrektan, #8 #9).

# CHANGES IN xfun VERSION 0.3

## NEW FEATURES

- Added a new functions `download_file()` to try various methods to download a file.

- Added a new function `is_ascii()` to test if a character vector only consists of ASCII characters.

- Added a new function `numbers_to_words()` to convert numbers to English words (thanks, @daijiang, #3).

# CHANGES IN xfun VERSION 0.2

## NEW FEATURES

- Added a `new_session` argument to `loadable()`.

- Added new functions `gsub_file()`, `gsub_files()`, `gsub_dir()`, and `gsub_ext()` to replace strings in files.

- Added new functions `Rscript` and `Rcmd` as wrappers of `system2('Rscript')` and `system2('R', 'CMD')`, respectively.

- Added a new function `install_dir()` to install a source package from a directory.

- Added a new function `file_string()` to read a text file (encoded in UTF-8) and return its content a single character string (lines concatenated by `\n`). 

- Added a new function `raw_string()` to print a character vector in its "raw" form using `cat(..., sep = '\n')` instead of `print()`, because the latter may introduce `[1]`, "extra" double quotes, and escape sequences, which are not very human-readable.

- Added a new function `session_info()` as an alternative to `sessionInfo()`.

- Added a new function `rev_check()` to run `R CMD check` on the reverse dependencies of a package, and a corresponding helper function `compare_Rcheck()` for showing the differences in logs with the CRAN version and the current version of the package, respectively.

- Added new functions for dealing with Markdown text: `prose_index()` returns the line indices of text that is prose (not code blocks), and `protect_math()` protects math expressions in Markdown in backticks.

- Added an `error` argument to `read_utf8()` to signal an error if the file is not encoded in UTF-8.

# CHANGES IN xfun VERSION 0.1

## NEW FEATURES

- `attr()` as an abbreviation of `base::attr(exact = TRUE)`.

- `file_ext()`, `sans_ext()`, and `with_ext()` to manipulate extensions in filenames.

- `in_dir()` to evaluate an R expression in a directory.

- `isFALSE()` as an abbreviation of `identical(x, FALSE)`.

- `is_windows()`, `is_macos()`, `is_linux()`, and `is_unix()` to test operating systems.

- `native_encode()` to try to encode a character vector in the native encoding.

- `normalize_path()` as an abbreviation of `normalizePath(winslash = '/', mustWork = FALSE)`.

- `optipng()` to run the command `optipng` to optimize all PNG files under a directory.

- `parse_only()` parses R code without keeping the source references.

- `pkg_attach()` and `pkg_load()` to attach and load a vector of packages, respectively (and optionally, install the missing packages).

- `read_utf8()` and `write_utf8()` to read and write UTF-8 files, respectively.

- `same_path()` to test if two paths are the same.

- `strict_list()` is a version of `list()` that disables partial matching of the `$` operator.

- `tojson()` is a simple JSON serializer.

- `try_silent()` is an abbreviation of `try(silent = TRUE)`.
