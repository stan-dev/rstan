
# cli 1.1.0

* cli has now functions to add ANSI styles to text. These use the crayon
  package internally, and provide a simpler interface. See the `col_*`,
  `bg_*`, `style_*` and also the `make_ansi_style()` and
  `combine_ansi_styles()` functions (#51).

* New `is_dynamic_tty()` function detects if `\r` should be used for a
  stream (#62).

* New `is_ansi_tty()` function detects if ANSI control sequences can be
  used for a stream.

* New `ansi_hide_cursor()`, `ansi_show_cursor()` and
  `ansi_with_hidden_cursor()` functions to hide and show the cursor in
  terminals.

* New `make_spinner()` function helps integrating spinners into your
  functions.

* Now `symbol` always uses ASCII symbols when the `cli.unicode` option is
  set to `FALSE`.

# 1.0.1

* New `cli_sitrep()` function, situation report about UTF-8 and ANSI
  color support (#53).

* Fall back to ASCII only characters on non-Windows platforms without
  UTF-8 support, and also in LaTeX when running knitr (#34).

# cli 1.0.0

First public release.
