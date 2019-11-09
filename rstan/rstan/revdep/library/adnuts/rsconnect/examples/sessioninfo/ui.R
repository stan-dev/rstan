library(shiny)

shinyUI(fluidPage(
  h3("URL"),
  verbatimTextOutput("urlInfo"),
  
  h3("Session"),
  verbatimTextOutput("sessionInfo"),
  
  h3("User"),
  verbatimTextOutput("userInfo")
))
