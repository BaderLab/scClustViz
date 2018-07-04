library(shiny)
library(RColorBrewer)
load("test.RData")

ui <- fixedPage(
  fixedRow(
    column(6,plotOutput("tsne",brush="tsneBrush")),
    column(6,
           actionButton("addCellsA","Set A: Add Cells",icon("plus")),
           actionButton("removeCellsA","Set A: Remove Cells",icon("minus")),
           hr(),
           actionButton("addCellsB","Set B: Add Cells",icon("plus")),
           actionButton("removeCellsB","Set B: Remove Cells",icon("minus"))
    )
  )
)

server <- function(input,output,session) {
  selectedSets <- reactiveValues(a=NULL,b=NULL)

  plot_tsne <- function() {
    par(mar=c(3,3,1,1),mgp=2:0)
    plot(dr_viz)
    points(dr_viz[selectedSets$a,],pch=19,col=brewer.pal(3,"PRGn")[1])
    points(dr_viz[selectedSets$b,],pch=19,col=brewer.pal(3,"PRGn")[3])
    points(dr_viz[rownames(dr_viz) %in% selectedSets$a & rownames(dr_viz) %in% selectedSets$b,],pch=19,col="red")
  }
  output$tsne <- renderPlot({ print(plot_tsne()) })
  
  
  currSel <- reactive(rownames(brushedPoints(as.data.frame(dr_viz),input$tsneBrush,xvar="tSNE_1",yvar="tSNE_2")))
  observeEvent(input$addCellsA,{ 
    selectedSets$a <- append(selectedSets$a,currSel()[!currSel() %in% selectedSets$a]) 
  })
  observeEvent(input$removeCellsA,{ 
    selectedSets$a <- selectedSets$a[!selectedSets$a %in% currSel()]
  })
  observeEvent(input$addCellsB,{ 
    selectedSets$b <- append(selectedSets$b,currSel()[!currSel() %in% selectedSets$b]) 
  })
  observeEvent(input$removeCellsB,{ 
    selectedSets$b <- selectedSets$b[!selectedSets$b %in% currSel()]
  })
  
  
  
  output$brushedCells <- renderPrint(paste0("|",currSel(),"|"))
  output$selected <- renderPrint(paste0("|",selectedSets[["a"]],"|"))
  
}

shinyApp(ui=ui,server=server)