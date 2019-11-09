# pkgbuild 1.0.3

* Tests which wrote to the package library are now skipped on CRAN.

* `build()` can now build a tar.gz file directly (#55)

# pkgbuild 1.0.2

* `build()` and `compile_dll()` gain a `register_routines` argument, to
  automatically register C routines with
  `tools::package_native_routines_registration_skeleton()` (#50)

* `build()` will now warn if trying to build packages on R versions <= 3.4.2 on
  Windows with a space in the R installation directory (#49)

* `build()` will now message if a build contains long paths, which are unsupported on windows
  (#48)

* `compile_dll()` no longer doubles output, a regression caused by the styling callback.
  (https://github.com/r-lib/devtools/issues/1877)

* `build()` output is now styled like that in the rcmdcheck package
  (https://github.com/r-lib/devtools/issues/1874).

* `build()` no longer sets compile flags (#46)

# pkgbuild 1.0.1

* `compile_dll()` now does not supply compiler flags if there is an existing
  user defined Makevars file.

* `local_build_tools()` function added to provide a deferred equivalent to
  `with_build_tools()`. So you can add rtools to the PATH until the end of a
  function body.

# pkgbuild 1.0.0

* Add metadata to support Rtools 3.5 (#38).

* `build()` only uses the `--no-resave-data` argument in `R CMD build`
  if the `--resave-data` argument wasn't supplied by the user
  (@theGreatWhiteShark, #26)

* `build()` now cleans existing vignette files in `inst/doc` if they exist. (#10)

* `clean_dll()` also deletes `symbols.rds` which is created when `compile_dll()`
  is run inside of `R CMD check`.

* First argument of all functions is now `path` rather than `pkg`.



