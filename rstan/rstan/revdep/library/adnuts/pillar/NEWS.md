# pillar 1.3.1

## Bug fixes

- Fix off-by-one error in distribution of empty space (#141).

## Visible changes

- `NA` in names is no longer escaped with backticks.
- Don't add quotes for pillars formatted with their `format()` method (tidyverse/tibble#448).

## Internal changes

- Update base type abbrevs to rlang 0.3.0 (#140, @lionel-).
- Tests work again in a 256-color terminal (#129).


# pillar 1.3.0

## Visible changes

- Unknown data types are formatted using `format()`, not `as.character()` (#120).

- Multi-tier colonnades can always fill the last tier, even if the width isn't a proper multiple of `getOption("width")`. (Example: `options(width = 80, tibble.width = 200)` will print a wide tibble in three tiers, each 80 characters wide, with a total width of 240 characters.)

- Fixed mixed formatting (showing some pillars with maximum, and some with minimum width). If a pillar's minimum width is smaller than `getOption("width")`, it is shown nevertheless, abbreviated with dots if necessary.

## Interface changes

- `format_type_sum()` gains `width` argument (#73).

## Performance improvements

- Printing large multi-tier colonnades is much faster, the code that distributes pillars over tiers uses a much simpler and much faster algorithm (tidyverse/tibble#422).

- Printing is now faster overall, because less work is done for formatting in "subtle" style (gray of a fixed level), and because `fansi::strip_sgr()` is used instead of `crayon::strip_style()`.

- Slightly faster printing of colonnades by reusing an intermediate result.

## Internal

- `pillar()` no longer adds backticks if `title` is non-syntactic.

- `colonnade()` supports data frames and matrices. When printing, each sub-column is shown individually, using a title that resembles the syntax used to access it. Also supports recursively nested data frames (with data frame or matrix columns).

- Added fuzz tests for character colonnades of varying widths.

- Use `fansi::substr_ctl()` in favor of `crayon::col_substr()`.


# pillar 1.2.3

- Eliminate CRAN check warning about undeclared withr dependency.
- More defensive test to address CRAN check failures on Solaris.
- `colonnade()` now handles pillars named `"sep"` (#115).
- `pillar_shaft.character()` gains `min_width` argument.


# pillar 1.2.2

- Whole numbers are printed without a decimal dot again. Numbers that are the result of a whole number divided by a power of 10 (subject to a tolerance to account for floating-point imprecision) are shown without trailing decimal zeros, even if these zeros are significant according to the `pillar.sigfig` option (#105).
- New `new_pillar_title()` and `new_pillar_type()` to support consistent output in `glimpse()` (#31).
- New `format_type_sum()` generic that allows overriding the formatting of the type summary in the capital (#73).
- The `digits.secs` option is respected when computing the width for date-time values (#102).


# pillar 1.2.1

Display
-------

- Turned off using subtle style for digits that are considered insignificant.  Negative numbers are shown all red.  Set the new option `pillar.subtle_num` to `TRUE` to turn it on again (default: `FALSE`).
- The negation sign is printed next to the number again (#91).
- Scientific notation uses regular digits again for exponents (#90).
- Groups of three digits are now underlined, starting with the fourth before/after the decimal point. This gives a better idea of the order of magnitude of the numbers (#78).
- Logical columns are displayed as `TRUE` and `FALSE` again (#95).
- The decimal dot is now always printed for numbers of type `numeric`. Trailing zeros are not shown anymore if all displayed numbers are whole numbers (#62).
- Decimal values longer than 13 characters always print in scientific notation.

Bug fixes
---------

- Numeric values with a `"class"` attribute (e.g., `Duration` from lubridate) are now formatted using `format()` if the `pillar_shaft()` method is not implemented for that class (#88).
- Very small numbers (like `1e-310`) are now printed corectly (tidyverse/tibble#377).
- Fix representation of right-hand side for `getOption("pillar.sigfig") >= 6` (tidyverse/tibble#380).
- Fix computation of significant figures for numbers with absolute value >= 1 (#98).

New functions
-------------

- New styling helper `style_subtle_num()`, formatting depends on the `pillar.subtle_num` option.


# pillar 1.1.0

- `NA` values are now shown in plain red, without changing the background color (#70).
- New options to control the output, with defaults that match the current behavior unless stated otherwise:
    - `pillar.sigfig` to control the number of significant digits, for highlighting and truncation (#72),
    - `pillar.subtle` to specify if insignificant digits should be printed in gray (#72),
    - `pillar.neg` to specify if negative digits should be printed in red,
    - `pillar.bold` to specify if column headers should be printed in bold (default: `FALSE`, #76),
    - `pillar.min_title_chars` to specify the minimum number of characters to display for each column name (default: 15 characters, #75).
- Shortened abbreviations for types: complex: cplx -> cpl, function: fun -> fn, factor: fctr -> fct (#71).
- Date columns now show sub-seconds if the `digits.secs` option is set (#74).
- Very wide tibbles now print faster (#85).


# pillar 1.0.1

- Work around failing CRAN tests on Windows.


# pillar 1.0.0

Initial release.

## User functions

    pillar(x, title = NULL, width = NULL, ...)
    colonnade(x, has_row_id = TRUE, width = NULL, ...)
    squeeze(x, width = NULL, ...)

## Functions for implementers of data types

    new_pillar_shaft_simple(formatted, ..., width = NULL, align = "left", min_width = NULL, na_indent = 0L)
    new_pillar_shaft(x, ..., width, min_width = width, subclass)
    new_ornament(x, width = NULL, align = NULL)
    get_extent(x)
    get_max_extent(x)

## Utilities

    dim_desc(x)
    style_na(x)
    style_neg(x)
    style_num(x, negative, significant = rep_along(x, TRUE))
    style_subtle(x)

## Testing helper

    expect_known_display(object, file, ..., width = 80L, crayon = TRUE)

## Own S3 methods

    pillar_shaft(x, ...) # AsIs, Date, POSIXt, character, default, list, logical, numeric
    type_sum(x) # AsIs, Date, POSIXct, data.frame, default, difftime, factor, ordered
    is_vector_s3(x) # Date, POSIXct, data.frame, default, difftime, factor, ordered
    obj_sum(x) # AsIs, POSIXlt, default, list
    extra_cols(x, ...) # squeezed_colonnade
