########## UI ########## 
ui <- fixedPage(
  fixedRow(titlePanel(paste("scClustViz -",dataTitle))),
  hr(),
  
  ######## Cluster Resolution Selection ########
  fixedRow(titlePanel("Cluster Resolution Selection"),
           column(6,
                  fixedRow(column(6,uiOutput("resSelect"),align="left"),
                           column(6,align="right",
                                  actionButton("go","View clusters at this resolution",icon("play")),
                                  actionButton("save","Save this resolution as default",icon("bookmark")))),
                  radioButtons("deType",NULL,list("# of marker genes per cluster"="deMarker",
                                                  "# of DE genes to nearest neighbouring cluster"="deNeighb"),inline=T),
                  plotOutput("cqPlot",height="500px")),
           column(6,plotOutput("sil",height="600px"))
  ),
  fixedRow(
    column(6,downloadButton("cqPlotSave","Save as PDF"),align="left"),
    column(6,downloadButton("silSave","Save as PDF"),align="right")
  ),
  hr(),
  
  ######## Cell-type Clusters ########
  fixedRow(titlePanel("Cell-type Clusters")),
  fixedRow(
    column(6,
           if (length(cellMarkers) > 0) {
             radioButtons("tsneLabels","Labels:",inline=T,
                          choices=list("Cluster numbers"="cn","Cluster annotations"="ca"))
           } else {
             radioButtons("tsneLabels","Labels:",inline=T,
                          choices=list("Cluster numbers"="cn"))
           },
           strong("Click point on plot below to select cluster")),
    column(3,selectInput("mdScatterX","x axis:",
                         choices=colnames(md)[!sapply(md,function(X) is.factor(X) | is.character(X))],
                         selected="total_counts"),align="left"),
    column(3,selectInput("mdScatterY","y axis:",
                         choices=colnames(md)[!sapply(md,function(X) is.factor(X) | is.character(X))],
                         selected="total_features"),align="left")
  ),
  fixedRow(
    column(6,plotOutput("tsne",height="570px",click="tsneClick")),
    column(6,plotOutput("mdScatter",height="570px"))
  ),
  fixedRow(
    column(6,align="left",downloadButton("tsneSave","Save as PDF")),
    column(6,align="right",downloadButton("mdScatterSave","Save as PDF"))
  ),
  hr(),
  
  fixedRow(
    column(6,selectInput("tsneMDcol","Metadata:",choices=colnames(md),
                         selected=grep("phase",colnames(md),value=T,ignore.case=T)[1])),
    column(3,selectInput("mdFactorData","Metadata (factor):",
                         choices=colnames(md)[sapply(md,function(X) is.factor(X) | is.character(X))],
                         selected=grep("phase",
                                       colnames(md)[sapply(md,function(X) is.factor(X) | is.character(X))],
                                       value=T,ignore.case=T)[1])),
    column(3,radioButtons("mdFactorRA","Factor counts per cluster:",inline=T,
                          choices=list("Absolute"="absolute","Relative"="relative")))
  ),
  fixedRow(
    column(6,plotOutput("tsneMD",height="570px")),
    column(6,plotOutput("mdFactor",height="570px"))
  ),
  fixedRow(
    column(6,align="left",downloadButton("tsneMDSave","Save as PDF")),
    column(6,align="right",downloadButton("mdFactorSave","Save as PDF"))
  ),
  hr(),
  
  ######## Cluster-wise Gene Stats #########
  fixedRow(titlePanel("Cluster-wise Gene Stats")),
  fixedRow(
    column(2,radioButtons("heatG","Heapmap Genes:",
                          choices=list("DE vs tissue average"="deTissue",
                                       "Marker genes"="deMarker",
                                       "DE vs neighbour"="deNeighb"))),
    column(2,uiOutput("DEclustSelect")),
    column(2,downloadButton("deGeneSave","Download gene list"),
           downloadButton("heatmapSave","Save as PDF"),align="right"), 
    column(6,uiOutput("DEgeneSlider"))
  ),
  fixedRow(plotOutput("heatmap",height="600px")),
  hr(),
  
  fixedRow(
    column(1,uiOutput("genePlotClustSelect")),
    column(6,if (length(cellMarkers) > 0) {
      radioButtons("cgLegend",inline=T,label="Highlighted genes:",
                   choices=c("Cell-type markers"="markers",
                             "Gene symbols (regex)"="regex",
                             "Top DE genes (from heatmap)"="heatmap"))
    } else {
      radioButtons("cgLegend",inline=T,label="Highlighted genes:",
                   choices=c("Gene symbols (regex)"="regex",
                             "Top DE genes (from heatmap)"="heatmap"))
    }),
    column(4,align="right",textInput("GOI","Gene symbols (regex)",demoRegex)),
    column(1,actionButton("GOIgo","Search",icon=icon("search")))
  ),tags$style(type='text/css', "button#GOIgo { margin-top: 25px; }"),
  fixedRow(
    plotOutput("clusterGenes",height="600px",click="cgClick"),
    downloadButton("clusterGenesSave","Save as PDF")
  ),
  hr(),
  
  fixedRow(
    column(4,uiOutput("cgSelect")),
    column(8,
           radioButtons("boxplotGene",inline=T,label="Gene of interest:",
                        choices=c("Click from plot above"="click",
                                  "From gene symbols (regex entry)"="regex")),
           checkboxGroupInput("bxpOpts",label=NULL,selected=c("sct","rnk"),inline=T,
                              choices=list("Include scatterplot"="sct",
                                           "Include gene rank"="rnk")))
  ),
  fixedRow(plotOutput("geneTest",height="500px"),
           downloadButton("geneTestSave","Save as PDF")
  ),
  hr(),
  
  ######## Distribution of genes of interest #########
  fixedRow(titlePanel("Distribution of Genes of Interest")),
  fixedRow(
    column(4,
           fixedRow(
             radioButtons("plotClust1",inline=T,label="Plot:",selected="goi",
                          choices=list("clusters"="clust","gene expression overlay"="goi")),
             checkboxInput("plotLabel1",label="Include cluster labels",value=T)
           ),
           fixedRow(
             column(9,textInput("GOI1",label="Gene symbols (regex):",demoRegex)),
             column(3,actionButton("GOI1go","Search",icon=icon("search")))
           ),tags$style(type='text/css', "button#GOI1go { margin-top: 25px; margin-left: -25px; }")
    ),
    column(2,uiOutput("GOI1select")),
    column(4,
           fixedRow(
             radioButtons("plotClust2",inline=T,label="Plot:",selected="goi",
                          choices=list("clusters"="clust","gene expression overlay"="goi")),
             checkboxInput("plotLabel2",label="Include cluster labels",value=T)
           ),
           fixedRow(
             column(9,textInput("GOI2",label="Gene symbols (regex):",demoRegex)),
             column(3,actionButton("GOI2go","Search",icon=icon("search")))
           ),tags$style(type='text/css', "button#GOI2go { margin-top: 25px; margin-left: -25px; }")
    ),
    column(2,uiOutput("GOI2select"))
  ),
  fixedRow(
    column(6,strong("If multiple genes are selected, the max expression per cell will be displayed")),
    column(6,strong("If multiple genes are selected, the max expression per cell will be displayed"))
  ),
  fixedRow(
    column(6,plotOutput("goiPlot1",height="600px")),
    column(6,plotOutput("goiPlot2",height="600px"))
  ),
  fixedRow(
    column(6,align="left",downloadButton("goiPlot1Save","Save as PDF")),
    column(6,align="right",downloadButton("goiPlot2Save","Save as PDF"))
  ),
  h1()
)

