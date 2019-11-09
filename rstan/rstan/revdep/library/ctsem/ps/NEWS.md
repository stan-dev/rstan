
# ps 1.3.0

* New `ps_cpu_count()` function returns the number of logical or
  physical processors.

# ps 1.2.1

* Fix a crash on Linux, that happened at load time (#50).

# ps 1.2.0

* New `ps_connections()` to list network connections. The
  `CleanupReporter()` testthat reporter can check for leftover open
  network connections in test cases.

* `ps_open_files()` does not include open sockets now on Linux, they are
  rather included in `ps_connections()`.

* `CleanupReporter()` now ignores `/dev/urandom`, some packages (curl,
  openssl, etc.) keep this file open.

* Fix `ps()` printing without the tibble package (#43).

* Fix compilation with ICC (#39).

* Fix a crash on Linux (#47).

# ps 1.1.0

* New `ps_num_fds()` returns the number of open files/handles.

* New `ps_open_files()` lists all open files of a process.

* New `ps_interrupt()` interrupts a process. It sends a `SIGINT` signal on
  POSIX systems, and it can send CTRL+C or CTRL+BREAK events on Windows.

* New `ps_users()` lists users connected to the system.

* New `ps_mark_tree()`, `ps_find_tree()`, `ps_kill_tree()`,
  `with_process_cleanup()`: functions to mark and clean up child
  processes.

* New `CleanupReporter`, to be used with testthat: it checks for
  leftover child processes and open files in `test_that()` blocks.

# ps 1.0.0

First released version.
