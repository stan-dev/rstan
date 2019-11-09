library(shiny)

shinyServer(function(input, output, session) {

  # output sessionInfo
  output$sessionInfo <- renderPrint({
    sessionInfo()
  })

  # output urlInfo
  output$urlInfo <- renderText({
    paste(sep = "",
          "protocol: ", session$clientData$url_protocol, "\n",
          "hostname: ", session$clientData$url_hostname, "\n",
          "pathname: ", session$clientData$url_pathname, "\n",
          "port: ",     session$clientData$url_port,     "\n",
          "search: ",   session$clientData$url_search,   "\n"
    )
  })

  # output userInfo
  output$userInfo <- renderText({
    paste(sep="",
          "user: ", session$user, "\n")
  })

})
