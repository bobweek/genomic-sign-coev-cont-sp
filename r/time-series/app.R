library(shiny)

ui <- shinyUI(fluidPage(
  mainPanel(
    htmlOutput("inc")
  )
))

#  ----- server.R -----

server <- function(input, output) {
  getPage<-function() {
    return(includeHTML("time-series-report.html"))
  }
  output$inc<-renderUI({getPage()})
}

shinyApp(ui, server)