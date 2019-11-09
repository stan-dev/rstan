## from Motoki Wu
library(R2admb)
setup_admb()

file.copy(file.path(Sys.getenv("ADMB_HOME"),
                    "examples","admb","simple","simple.tpl"),"./simple.tpl")
x <- rnorm(20, 5, 2)
y <- rnorm(20, 2, 9)
(d <- do_admb("simple",
              data = list(nobs = length(x), Y = y, x = x),
              params = list(a = 2, b = 1),
              run.opts=run.control(checkparam = "ignore", checkdata = "ignore")))
summary(d)
unlink(c("simple","simple.tpl"))
