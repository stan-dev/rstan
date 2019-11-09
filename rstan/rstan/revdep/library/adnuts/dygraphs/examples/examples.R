
library(dygraphs)

hw <- HoltWinters(ldeaths)
p <- predict(hw, n.ahead = 36, prediction.interval = TRUE, level = 0.95)
ts <- cbind(ldeaths, p)

dygraph(ts, "Deaths from Lung Disease (UK)") %>%
  dySeries("ldeaths", label = "Deaths") %>%
  dySeries(c("p.lwr", "p.fit", "p.upr"), label = "Predicted") %>%
  dyRangeSelector()

