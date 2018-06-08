ui <- pageWithSidebar(
  headerPanel('Iris k-means clustering'),
  sidebarPanel(
    uiOutput("test1"),
    uiOutput("test2")
  ),
  mainPanel(
    textOutput("selected")
  )
)

server <- function(input, output, session) {
  selected <- reactiveValues(cl="")

  output$test1 <- renderUI({
    selectInput("test1","Cluster:",choices=c("",1:3),selected=selected$cl)
  })
  output$test2 <- renderUI({
    selectInput("test2","Cluster:",choices=c("",1:3),selected=selected$cl)
  })
  
  observeEvent(input$test1,{selected$cl <- input$test1})
  observeEvent(input$test2,{selected$cl <- input$test2})
  
  output$selected<- renderText(paste0("|",selected$cl,"|"))
  
}

shinyApp(ui=ui,server=server)
