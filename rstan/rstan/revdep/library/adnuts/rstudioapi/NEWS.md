# rstudioapi 0.10

* Added the parameters `echo` and `focus` to `sendToConsole()`.

# rstudioapi 0.9

* Added functions for displaying jobs in RStudio's Jobs pane: `jobAdd()`, `jobRemove()`, etc.

* Added `translateLocalUrl()`, for translating localhost URLs to externally addressable ones on RStudio Server.

# rstudioapi 0.8

* Added functions for installing + using build tools:
  `buildToolsCheck()`, `buildToolsInstall()`, `buildToolsExec()`
  
* Added functions for installing + using themes: `addTheme()`, `applyTheme()`,
  `convertTheme()`, `getThemes()`, `getThemeInfo()`.

* Added `previewSql()`, for previewing output from executing a SQL query.

* Added `askForSecret()`, for prompting the user to enter a password or otherwise privileged information.

* Fixed an issue where `getActiveProject()` failed for non-ASCII paths. (#86)

# rstudioapi 0.7

* Added methods for prompting the user for file paths: `selectFile()`,
  `selectDirectory()`.

* `askForPassword()` gains a default prompt (#41)

* Add `createProjectTemplate()` function

* Add `setPersistentValue()` / `getPersistentValue()` functions

* Add methods for interacting with Terminal tab:
  `terminalActivate()`, `terminalClear()`, `terminalCreate()`, `terminalList()`,
  `terminalBuffer()`, `terminalContext()`, `terminalVisible()`, `terminalBusy()`,
  `terminalRunning()`, `terminalKill()`, `terminalSend()`, `terminalExecute()`,
  and `terminalExitCode()`.

# rstudioapi 0.6

* Add sendToConsole function

* Add APIs for setting cursor position in document

# rstudioapi 0.5

* Add askForPassword function

* Add getActiveProject function

# rstudioapi 0.4

* Add API methods for interacting with a document open in RStudio: 'insertText()', 'modifyRange()' and 'getActiveDocumentContext()'.

# rstudioapi 0.3

* Add stub and documentation for sourceMarker function

# rstudioapi 0.2

* Compatibility with calling conventions for RStudio v0.99

* Stubs and documentation for versionInfo, previewRd, and viewer functions

# rstudioapi 0.1

Initial release to CRAN

