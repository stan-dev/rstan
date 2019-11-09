
## -------------------------------------------------------------

.l

app <- cliapp$new()

f1 <- function() {
  app$h1("Title")
  app$text(lorem_ipsum())
  
  app$h2("Subtitle")
  app$text(lorem_ipsum())
  
  app$h3("Subsubtitle")
  app$verbatim("foo\n", "bar\n", "foobar\n")
}

## -------------------------------------------------------------

f2 <- function() {
  app$alert_success("All is good")
  app$alert_info("For your information")
  app$alert_warning("Something might be wrong")
  app$alert_danger("This is definitely wrong")
}

## -------------------------------------------------------------

f3 <- function() {
  x <- 1
  app$h1("Header 1")
  app$text("glue substitution is automatic: {x}")
  app$h1("{x} Even in Headers")
}

## -------------------------------------------------------------

f4 <- function() {
  app$h1("About the {pkg callr} package")
  app$text("code: {code sapply(x, f)} function: {fun callr::r}")
}

## -------------------------------------------------------------

f5 <- function() {  
  app$par()
  app$text(lorem_ipsum())
  
  app$par()
  app$text(lorem_ipsum())
}

## -------------------------------------------------------------

f6 <- function() {
  app$ul()
  app$it("foo")
  app$it(c("bar", "foobar"))

  lid <- app$ul()
  app$it("sublist")
  app$it(c("sub", "subsub"))

  app$ol()
  app$it("First")
  app$it(c("Second", "Third"))

  app$end(lid)
  app$it("bar again!")
}

## f5()
## f5()

## -------------------------------------------------------------

## Demo CSS styling

red <- list(
  ".red" = list(
    color = "red",
    "background-color" = "yellow"
  )
)
options(cli.theme = red)

app <- cliapp$new()

f7 <- function() {
  app$par()
  app$text(lorem_ipsum())

  app$par(class = "red")
  app$text(lorem_ipsum())
  app$end()
  
  app$par()
  app$text(lorem_ipsum())
}

## -------------------------------------------------------------

f8 <- function() {
  x <- 100
  app$alert_success("so far so good: {x}")
  bar <- app$progress_bar(total = 5)
  bar$tick()
  Sys.sleep(1/2)
  bar$tick()
  Sys.sleep(1/2)

  app$alert_success("still very good: {x}!")
  Sys.sleep(1)
  bar$tick()
  Sys.sleep(1/2)
  app$text(lorem_ipsum())
  bar$tick()
  Sys.sleep(1)
  bar$tick()
  Sys.sleep(1/2)
  app$alert_success("aaaaand we are done")
}
