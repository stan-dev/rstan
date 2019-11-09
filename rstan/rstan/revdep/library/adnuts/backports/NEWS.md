# backports 1.1.4
* Fixed import of `warningCondition` and `errorCondition`.

# backports 1.1.3

* Added `warningCondition()` and `errorCondition()` for R versions prior to 3.6.0.
* Added `capture.output()` with support for argument `type` for R versions prior to 3.3.0.
* Added `URLencode()` with support for argument `repeated` for R versions prior to 3.2.0.

# backports 1.1.2

* Improved import mechanism.
* Added `.valid.factor()` for R versions prior to 3.4.0.

# backports 1.1.1

* Added `...length()` and `...elt()` for R versions prior to 3.5.0.
* Added `isFALSE()` for R versions prior to 3.5.0.

# backports 1.1.0

* New import mechanism to import packages during load-time with the function `import()`.
  This is now the recommended way to use backports non-interactively.
  Simply importing backports in the NAMESPACE still works, but is comparably error-prone
  if the same library is used by multiple R installations.

# backports 1.0.5

* Added `get0()` for R versions prior to 3.2.0.
* Added examples.

# backports 1.0.4

* Added `hasName()` for R versions prior to 3.4.0.
* Added `file.info()` with backport for argument `extra_cols`.

# backports 1.0.3

* Removed stringi dependency.

# backports 1.0.2

* Fixed `file.size()`, `file.mtime()` and `file.mode()` for R-3.1.x.

# backports 1.0.1

* Added `file.size()`, `file.mtime()` and `file.mode()` for R versions prior to 3.2.0.

# backports 1.0.0

* Initial version.
