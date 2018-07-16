########## UI ########## 
ui <- fixedPage(
  fixedRow(
    titlePanel(paste("scClustViz -",dataTitle)),
    includeMarkdown(paste0(dataPath,"intro.md"))
    ),
  hr(),
  
  ######## Clustering Solution Selection ########
  fixedRow(
    titlePanel("Clustering Solution Selection"),
    p(paste("Here you can compare the results of clustering at different resolutions to",
            "determine the appropriate clustering solution for your data. You can see the",
            "cluster solutions represented as boxplots on the left, where each boxplot",
            "represents the number of genes differentially expressed between each cluster",
            "and its nearest neighbour, or marker genes per cluster. The cluster selected",
            "in the pulldown menu is highlighted in red, and the silhouette plot for that",
            "cluster is shown on the right.")),
    p(paste("A silhouette plot is a horizontal barplot where each bar is a cell, grouped by",
            "cluster. The width of each bar represents the difference between mean distance",
            "to other cells within the cluster and mean distance to cells in the nearest",
            "neighbouring cluster. Distance is Euclidean in reduced dimensional space.",
            "Positive silhouettes indicate good cluster cohesion.")),
    p(paste("Once you've selected an appropriate cluster solution (we suggest picking one",
            "where all nearest neighbouring clusters have differentially expressed genes",
            "between them), click 'View clusters at this resolution' to proceed. If you",
            "want to save this cluster solution as the default for next time, click 'Save",
            "this resolution as default'. All figures can be downloaded as PDFs by clicking",
            "the buttons next to each figure.")),
    h1()
  ),
  fixedRow(
    column(6,
           fixedRow(column(6,uiOutput("resSelect"),align="left"),
                    column(6,align="right",
                           actionButton("go","View clusters at this resolution",icon("play")),
                           uiOutput("saveButton")
                    )
           ),
           radioButtons("deType",NULL,list("# of DE genes to nearest neighbouring cluster"="deNeighb",
                                           "# of marker genes per cluster"="deMarker"),inline=T),
           plotOutput("cqPlot",height="500px")),
    column(6,plotOutput("sil",height="600px"))
  ),
  fixedRow(
    column(6,downloadButton("cqPlotSave","Save as PDF"),align="left"),
    column(6,downloadButton("silSave","Save as PDF"),align="right")
  ),
  hr(),
  
  ######## Dataset and Cluster Metadata Inspection ########
  fixedRow(
    titlePanel("Dataset and Cluster Metadata Inspection"),
    p(paste("Here you can explore your dataset as a whole: cluster assignments for all",
            "cells; metadata overlays for cell projections; and figures for comparing",
            "both numeric and categorical metadata.")),
    p(paste("You can select any cluster for further assessment by clicking on a cell",
            "from that cluster in the first figure on the left.")),
    h1()
  ),
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
    column(6,selectInput("tsneMDcol","Metadata:",choices=colnames(md),
                         selected=grep("phase",colnames(md),value=T,ignore.case=T)[1]))
  ),
  fixedRow(
    column(6,plotOutput("tsne",height="580px",click="tsneClick")),
    column(6,plotOutput("tsneMD",height="580px"))
  ),
  fixedRow(
    column(6,align="left",downloadButton("tsneSave","Save as PDF")),
    column(6,align="right",downloadButton("tsneMDSave","Save as PDF"))
  ),
  hr(),
  
  fixedRow(
    column(3,selectInput(
      "mdScatterX","x axis:",
      choices=colnames(md)[!sapply(md,function(X) is.factor(X) | is.character(X))],
      selected="total_counts"),align="left"),
    column(3,selectInput(
      "mdScatterY","y axis:",
      choices=colnames(md)[!sapply(md,function(X) is.factor(X) | is.character(X))],
      selected="total_features"),align="left"),
    column(3,selectInput(
      "mdFactorData","Metadata (factor):",
      choices=colnames(md)[sapply(md,function(X) is.factor(X) | is.character(X))],
      selected=grep("phase",
                    colnames(md)[sapply(md,function(X) is.factor(X) | is.character(X))],
                    value=T,ignore.case=T)[1])),
    column(3,radioButtons("mdFactorRA","Factor counts per cluster:",inline=T,
                          choices=list("Absolute"="absolute","Relative"="relative")))
  ),
  fixedRow(
    column(6,plotOutput("mdScatter",height="560px")),
    column(6,plotOutput("mdFactor",height="560px"))
  ),
  fixedRow(
    column(6,align="left",downloadButton("mdScatterSave","Save as PDF")),
    column(6,align="right",downloadButton("mdFactorSave","Save as PDF"))
  ),
  hr(),
  
  ######## Cluster-wise Gene Stats #########
  fixedRow(
    titlePanel("Cluster-wise Gene Stats"),
    p(paste()),
    h1()
    
  ),
  fixedRow(
    column(2,uiOutput("heatDEtype")),
    column(2,uiOutput("DEclustSelect")),
    column(2,downloadButton("deGeneSave","Download gene list"),
           downloadButton("heatmapSave","Save as PDF"),align="right"), 
    column(6,uiOutput("DEgeneSlider"))
  ),
  fixedRow(plotOutput("heatmap",height="600px")),
  hr(),
  
  fixedRow(
    column(2,uiOutput("genePlotClustSelect")),
    column(6,if (length(cellMarkers) > 0) {
      radioButtons("cgLegend",inline=T,label="Highlighted genes:",
                   choices=c("Cell-type markers"="markers",
                             "Top DE genes (from heatmap)"="heatmap",
                             "Gene symbols (regex)"="regex"))
    } else {
      radioButtons("cgLegend",inline=T,label="Highlighted genes:",
                   choices=c("Top DE genes (from heatmap)"="heatmap",
                             "Gene symbols (regex)"="regex"))
    }),
    column(3,align="right",textInput("GOI","Gene symbols (regex)",demoRegex)),
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
  hr(),
  
  ######## Custom sets for DE #########
  fixedRow(titlePanel("Select cells for direct comparison")),
  fixedRow(
    column(6,
           plotOutput("tsneSelDE",brush="tsneBrush",height="580px")),
    column(6,
           p(paste("Here you can select cells to further explore using the figures above.",
                   "Click and drag to select cells, and use the buttons below to add them",
                   "to a set of cells. When your sets are ready, name the comparison and",
                   "click the 'Calculate & Save' button. Once the calculation is done",
                   "the comparison will be added to the cluster list at the top of the page.")),
           hr(),
           selectInput("tsneSelDEcol","Metadata overlay:",choices=c("",colnames(md))),
           hr(),
           column(6,htmlOutput("textSetA"),
                  actionButton("addCellsA","Set A: Add Cells",icon("plus"),
                               style="color: #fff; background-color: #a50026"),
                  actionButton("removeCellsA","Set A: Remove Cells",icon("minus"),
                               style="color: #a50026; background-color: #fff; border-color: #a50026")
           ),
           column(6,htmlOutput("textSetB"),
                  actionButton("addCellsB","Set B: Add Cells",icon("plus"),
                               style="color: #fff; background-color: #313695"),
                  actionButton("removeCellsB","Set B: Remove Cells",icon("minus"),
                               style="color: #313695; background-color: #fff; border-color: #313695")
           ),
           htmlOutput("textOverlap"),
           hr(),
           textInput("DEsetName","Short name for this comparison:",
                     placeholder="A-z0-9 only please"),
           actionButton("calcDE","Calculate differential gene expression",icon("play")),
           hr(),
           span(textOutput("calcText"),style="color:red")
    )
  ),
  h1()
)


########## Server ##########
server <- function(input,output,session) {
  d <- reactiveValues(cl=cl,CGS=CGS,
                      clusterID=clusterID,
                      deTissue=deTissue,
                      deMarker=deMarker)
  
  clustCols <- reactive({
    if (grepl("^Comp",input$res)) {
      c(brewer.pal(3,"PRGn")[c(1,3)],"grey80")
    } else if (length(levels(d$cl[,input$res])) <= 8) {
      brewer.pal(length(levels(d$cl[,input$res])),"Dark2")[1:length(levels(d$cl[,input$res]))]
    } else {
      rainbow2(length(levels(d$cl[,input$res])))
    }
  })
  
  
  ######## Cluster Resolution Selection ########
  #### Inter-cluster DE boxplots ####
  numClust <- sapply(cl[,!grepl("^Comp",colnames(cl))],function(X) length(levels(X)))
  clustList <- reactive({ 
    temp <- as.list(colnames(d$cl)) 
    names(temp)[seq_along(numClust)] <- paste0(unlist(temp)[seq_along(numClust)],
                                               ": ",numClust," clusters")
    if (length(temp) > length(numClust)) {
      names(temp)[setdiff(seq_along(temp),seq_along(numClust))] <- 
        paste0("Comparison: ",
               sub("Comp.","",fixed=T,
                   x=unlist(temp)[setdiff(seq_along(temp),seq_along(numClust))]))  
    }
    return(temp)
  })
  output$resSelect <- renderUI({
    selectInput("res","Resolution:",choices=clustList(),selected=savedRes)
  })
  output$saveButton <- renderUI({
    if (grepl("^Comp",input$res)) {
      actionButton("updateForViz","Save this comparison to disk",icon("save"))
    } else {
      actionButton("save","Save this resolution as default",icon("bookmark"))
    }
  })
  numClust <- numClust[numClust > 1]
  
  plot_cqPlot <- function() {
    numDEgenes <- lapply(get(input$deType)[!grepl("^Comp",names(get(input$deType)))],
                         function(X) sapply(X,nrow))
    toplim <- c(21,max(unlist(numDEgenes)) + 20)
    botlim <- c(-1,21)
    
    if (grepl("^Comp",input$res)) {
      par(mar=c(3,3.5,1,1))
      plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
      text(.5,.5,paste("Press 'View clusters at this resolution'",
                       "to view the comparison",
                       sub("Comp.","",input$res,fixed=T),sep="\n"))
    } else {
      par(mar=c(0.2,3.5,1,1),mgp=2:0,mfrow=c(2,1))
      plot(x=numClust,y=sapply(numDEgenes,median),type="l",
           xlim=range(numClust)+c(-.5,.5),ylim=toplim,yaxs="i",xaxt="n",ylab=NA)
      abline(h=seq(0,max(unlist(numDEgenes)),10),lty=3,col=alpha(1,0.3))
      for (i in names(numDEgenes)[names(numDEgenes) != input$res]) {
        boxplot(numDEgenes[[i]],add=T,at=numClust[i],yaxt="n")
      }
      if (any(names(numDEgenes) == input$res)) {
        boxplot(numDEgenes[[input$res]],add=T,at=numClust[input$res],border="red")
      }
      
      par(mar=c(3,3.5,0.2,1),mgp=2:0)
      plot(x=numClust,y=sapply(numDEgenes,median),type="l",
           xlim=range(numClust)+c(-.5,.5),ylim=botlim,yaxs="i",xlab="Number of clusters",ylab=NA)
      abline(h=seq(0,max(unlist(numDEgenes)),10),lty=3,col=alpha(1,0.3))
      for (i in names(numDEgenes)[names(numDEgenes) != input$res]) {
        boxplot(numDEgenes[[i]],add=T,at=numClust[i],yaxt="n")
      }
      if (any(names(numDEgenes) == input$res)) {
        boxplot(numDEgenes[[input$res]],add=T,at=numClust[input$res],border="red")
      }
      mtext(switch(input$deType,
                   "deMarker"="Positive DE genes per cluster to all other clusters",
                   "deNeighb"="Positive DE genes per cluster to nearest cluster")
            ,side=2,line=2.5,at=botlim[2],xpd=NA)
    }
  }
  
  output$cqPlot <- renderPlot({
    print(plot_cqPlot())
  })
  
  output$cqPlotSave <- downloadHandler(
    filename="cqPlot.pdf",
    content=function(file) {
      pdf(file,width=7,height=6)
      print(plot_cqPlot())
      dev.off()
    }
  )
  
  #### Silhouette plot ####
  plot_sil <- function() {
    tempSil <- silhouette(as.integer(d$cl[,input$res]),dist=silDist)
    par(mar=c(4.5,.5,1.5,1.5),mgp=2:0)
    if (length(tempSil) <= 1) {
      plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
      text(.5,.5,paste("Silhouette plot cannot be computed",
                       "with less than two clusters.",sep="\n"))
    } else {
      plot(tempSil,beside=T,border=NA,main=NA,col=clustCols(),do.n.k=T)
    }
  }
  
  output$sil <- renderPlot({
    print(plot_sil())
  })
  
  output$silSave <- downloadHandler(
    filename="sil.pdf",
    content=function(file) {
      pdf(file,width=6,height=7)
      print(plot_sil())
      dev.off()
    }
  )
  
  #### res buttons ####
  res <- eventReactive(input$go,input$res,ignoreNULL=F)
  
  observeEvent(input$save,{
    savedRes <<- input$res #<<- updates variable outside scope of function (ie. global environment)
    save(savedRes,file=paste0(dataPath,dataTitle,"_savedRes.RData"))
  })
  
  
  ######## Cell-type Clusters ########
  clusts <- reactive(d$cl[,res()])
  
  #### Cell-type tSNE ####
  plot_tsne_labels <- function() {
    if (input$tsneLabels == "ca") {
      temp_labelNames <- sapply(unique(d$clusterID[[res()]]),function(X) 
        names(which(d$clusterID[[res()]] == X)),simplify=F)
      temp_labels <- apply(dr_viz,2,function(Y) 
        tapply(Y,apply(sapply(temp_labelNames,function(X) clusts() %in% X),1,which),mean))
      if (!is.matrix(temp_labels)) { temp_labels <- rbind(temp_labels) }
      text(temp_labels,labels=names(temp_labelNames),font=2,cex=1.5)
    } else if (input$tsneLabels == "cn") {
      temp_labels <- apply(dr_viz,2,function(X) tapply(X,clusts(),mean))
      if (!is.matrix(temp_labels)) { temp_labels <- rbind(temp_labels) }
      text(temp_labels,labels=levels(clusts()),font=2,cex=1.5)
    } else {
      legend("center",legend="You changed the label choice names...")
    }
  }
  
  plot_tsne <- function() {
    par(mar=c(3,3,4,1),mgp=2:0)
    plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
         main=paste("tSNE at",res(),"using",ncol(dr_clust),"PCs"),
         xlim=range(dr_viz[,1]),ylim=range(dr_viz[,2]))
    if (any(ci())) {
      points(dr_viz[!ci(),],pch=21,
             col=alpha(clustCols()[clusts()],0.2)[!ci()],
             bg=alpha(clustCols()[clusts()],0.1)[!ci()])
      points(dr_viz[ci(),],pch=21,
             col=alpha(clustCols()[clusts()],1)[ci()],
             bg=alpha(clustCols()[clusts()],0.5)[ci()])
    } else {
      points(dr_viz,pch=21,
             col=alpha(clustCols()[clusts()],1),
             bg=alpha(clustCols()[clusts()],0.5))
    }
    if (hiC() != "") {
      mtext(side=3,line=-1,text=paste("Cluster",hiC(),"-",
                                      d$clusterID[[res()]][hiC()],"-",
                                      sum(clusts() == hiC()),"cells"))
    }
  }
  
  output$tsne <- renderPlot({
    if (length(res()) > 0) {
      print(plot_tsne())
      print(plot_tsne_labels())
    }
  })
  
  output$tsneSave <- downloadHandler(
    filename="tsne.pdf",
    content=function(file) {
      pdf(file,width=7,height=7)
      print(plot_tsne())
      print(plot_tsne_labels())
      dev.off()
    }
  )
  
  #### clusterSelect ####
  
  
  clusterSelect <- reactiveValues(cl=NULL)

  observeEvent(input$tsneClick,{ clusterSelect$cl <- input$tsneClick })

  cSelected <- reactive({
    t <- nearPoints(as.data.frame(dr_viz),clusterSelect$cl,xvar="tSNE_1",yvar="tSNE_2",threshold=5)
    t2 <- d$cl[rownames(t)[1],res()]
    if (is.na(t2)) { 
      return("") 
    } else if (t2 == "Unselected") {
      return("")
    } else { 
      return(t2) 
    }
  })
  
  hiC <- reactive({ 
    if (length(res()) < 1) {
      return("")
    } else if (input$genePlotClust != "") {
      d$cl[which(d$cl[,res()] == input$genePlotClust)[1],res()]
    } else {
      return(input$genePlotClust)
    }
  })
  
  ci <- reactive({
    if (hiC() == "") {
      rep(F,length(clusts()))
    } else {
      clusts() == hiC()
    }
  })
  
  #### Metadata tSNE overlay ####
  plot_tsneMD <- function() {
    if (is.factor(md[,input$tsneMDcol]) | is.character(md[,input$tsneMDcol])) {
      id <- as.factor(md[,input$tsneMDcol])
      if (length(levels(md[,input$tsneMDcol])) <= 8) {
        idcol <- brewer.pal(length(levels(md[,input$tsneMDcol])),
                            "Dark2")[1:length(levels(md[,input$tsneMDcol]))]
      } else {
        idcol <- rainbow2(length(levels(md[,input$tsneMDcol])))
      }
    } else {
      id <- cut(md[,input$tsneMDcol],100)
      idcol <- viridis(100,d=-1)
    }
    layout(cbind(2:1),heights=c(1,9))
    par(mar=c(3,3,0,1),mgp=2:0)
    plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
         xlim=range(dr_viz[,1]),ylim=range(dr_viz[,2]))
    if (any(ci())) {
      points(dr_viz[!ci(),],pch=21,
             col=alpha(idcol,.1)[id[!ci()]],
             bg=alpha(idcol,0.05)[id[!ci()]])
      points(dr_viz[ci(),],pch=21,
             col=alpha(idcol,.8)[id[ci()]],
             bg=alpha(idcol,0.4)[id[ci()]])
    } else {
      points(dr_viz,pch=21,
             col=alpha(idcol,.8)[id],
             bg=alpha(idcol,0.4)[id])
    }
    plot_tsne_labels()
    if (is.factor(md[,input$tsneMDcol]) | is.character(md[,input$tsneMDcol])) {
      par(mar=c(0,0,0,0))
      plot.new()
      legend("bottom",bty="n",horiz=T,pch=c(NA,rep(21,length(levels(md[,input$tsneMDcol])))),
             legend=c(paste0(input$tsneMDcol,":"),levels(md[,input$tsneMDcol])),
             col=c(NA,idcol),pt.bg=c(NA,alpha(idcol,0.5)))
    } else {
      par(mar=c(0,5,3,3))
      barplot(rep(1,100),space=0,col=idcol,xaxt="n",yaxt="n",border=NA,main=input$tsneMDcol)
      text(x=c(1,100),y=1,pos=c(2,4),xpd=NA,labels=round(range(md[,input$tsneMDcol]),2))
    }
  }
  
  output$tsneMD <- renderPlot({
    if (length(res()) > 0) {
      print(plot_tsneMD())
    }
  })
  
  output$tsneMDSave <- downloadHandler(
    filename="tsneMD.pdf",
    content=function(file) {
      pdf(file,width=7,height=7)
      print(plot_tsneMD())
      dev.off()
    }
  )
  
  #### Metadata Factor Barplot ####
  plot_mdFactor <- function() {
    id <- switch(input$mdFactorRA,
                 "relative"=tapply(md[,input$mdFactorData],clusts(),
                                   function(X) table(X) / length(X)),
                 "absolute"=tapply(md[,input$mdFactorData],clusts(),table))
    if (is.list(id)) { id <- do.call(cbind,id) }
    idylab <- switch(input$mdFactorRA,
                     "relative"="Proportion of cells per cluster",
                     "absolute"="Number of cells per cluster")
    if (length(levels(md[,input$mdFactorData])) <= 8) {
      idcol <- brewer.pal(length(levels(md[,input$mdFactorData])),
                          "Dark2")[1:length(levels(md[,input$mdFactorData]))]
    } else {
      idcol <- rainbow2(length(levels(md[,input$mdFactorData])))
    }
    par(mar=c(3,3,4,1),mgp=2:0)
    barplot(id,col=idcol,ylab=idylab,
            legend.text=levels(md[,input$mdFactorData]),
            args.legend=list(x="topright",horiz=T,inset=c(0,-.08),bty="n"))
    mtext(input$mdFactorData,side=3,adj=0,font=2,line=1,cex=1.2)
  }
  
  output$mdFactor <- renderPlot({
    if (length(res()) > 0) {
      print(plot_mdFactor())
    }
  })
  
  output$mdFactorSave <- downloadHandler(
    filename="mdFactor.pdf",
    content=function(file) {
      pdf(file,width=7,height=7)
      print(plot_mdFactor())
      dev.off()
    }
  )
  
  #### Metadata Scatterplot ####
  plot_mdScatter <- function() {
    layout(matrix(c(2,1,0,3),2),c(5,1),c(1,5))
    par(mar=c(3,3,0,0),mgp=2:0,cex=1.1)
    if (all(ci())) {
      plot(md[,input$mdScatterX],md[,input$mdScatterY],
           pch=21,col=alpha("red",0.4),bg=alpha("red",0.2),
           xlab=input$mdScatterX,ylab=input$mdScatterY)
    } else {
      plot(md[!ci(),input$mdScatterX],md[!ci(),input$mdScatterY],
           pch=21,col=alpha("black",0.2),bg=alpha("black",0.1),
           xlab=input$mdScatterX,ylab=input$mdScatterY)
      points(md[ci(),input$mdScatterX],md[ci(),input$mdScatterY],
             pch=21,col=alpha("red",0.4),bg=alpha("red",0.2))
    }
    if (any(ci())) {
      legend("topleft",bty="n",pch=21,col="red",pt.bg=alpha("red",0.5),
             legend=paste("Cluster",hiC(),"-",d$clusterID[[res()]][hiC()]))
    }
    par(mar=c(0,3,1,0))
    boxplot(tapply(md[,input$mdScatterX],ci(),c),
            horizontal=T,xaxt="n",yaxt="n",border=c("black","red"))
    par(mar=c(3,0,0,1))
    boxplot(tapply(md[,input$mdScatterY],ci(),c),
            horizontal=F,xaxt="n",yaxt="n",border=c("black","red"))
  }
  
  output$mdScatter <- renderPlot({
    if (length(res()) > 0) {
      print(plot_mdScatter())
    }
  })
  
  output$mdScatterSave <- downloadHandler(
    filename="mdScatter.pdf",
    content=function(file) {
      pdf(file,width=7,height=7)
      print(plot_mdScatter())
      dev.off()
    }
  )
  
  
  ######## Cluster-wise Gene Stats #########
  
  #### Heatmap genes ####
  output$heatDEtype <- renderUI({
    if (grepl("^Comp",input$res)) {
      temp <- list("DE vs rest"="deTissue",
                   "Set A vs Set B"="deMarker")
    } else {
      temp <- list("DE vs rest"="deTissue",
                   "Marker genes"="deMarker",
                   "DE vs neighbour"="deNeighb")
    }
    radioButtons("heatG","Heapmap Genes:",choices=temp)
  })
  
  output$DEgeneSlider <- renderUI({
    if (length(res()) > 0) {
      switch(
        input$heatG,
        deTissue=
          sliderInput("DEgeneCount",min=1,max=max(sapply(d$deTissue[[res()]],nrow)),
                      value=5,step=1,ticks=T,width="100%",
                      label=HTML(paste(
                        "Positive differential gene expression of cluster over tissue",
                        "# of genes per cluster to show",sep="<br/>"
                      ))),
        deMarker=
          sliderInput("DEgeneCount",min=1,max=max(sapply(d$deMarker[[res()]],nrow)),
                      value=5,step=1,ticks=T,width="100%",
                      label=HTML(paste(
                        "Positive differential gene expression between cluster and all other clusters",
                        "# of genes per cluster to show",sep="<br/>"
                      ))),
        deNeighb=
          sliderInput("DEgeneCount",min=1,max=max(sapply(deNeighb[[res()]],nrow)),
                      value=5,step=1,ticks=T,width="100%",
                      label=HTML(paste(
                        "Positive differential gene expression between cluster and nearest neighbour",
                        "# of genes per cluster to show",sep="<br/>"
                      ))))
    }
  })
  
  output$DEclustSelect <- renderUI({
    if (length(res()) > 0) {
      selectInput("DEclustNum","Cluster # for gene list",
                  choices=levels(clusts())[!levels(clusts()) == "Unselected"])
    }
  })
  
  heatGenes <- reactive({
    temp <- unique(unlist(lapply(
      switch(input$heatG,
             deTissue=d$deTissue[[res()]],
             deMarker=d$deMarker[[res()]],
             deNeighb=deNeighb[[res()]]),
      function(X) 
        if (nrow(X) == 0) { NA } else { rownames(X)[1:input$DEgeneCount] }
    )))
    temp <- temp[!is.na(temp)]
    return(temp)
  })
  
  clustMeans <- reactive({ #This only works if input is in ascending order of adjusted p value.
    temp <- sapply(d$CGS[[res()]],function(X) X[heatGenes(),"MTC"])
    rownames(temp) <- heatGenes()
    return(t(temp))
  })
  
  hC <- reactive(hclust(dist(clustMeans()),"single"))
  hG <- reactive(hclust(dist(t(clustMeans())),"complete"))
  
  sepClust <- reactive({
    if (hiC() == "") {
      return(c(NA,NA))
    } else {
      return(nrow(clustMeans()) - 
               c(which(levels(clusts())[hC()$order] == hiC()) - 1,
                 which(levels(clusts())[hC()$order] == hiC())))
    }
  })
  
  plot_heatmap <- function() {
    if (length(levels(clusts())) <= 1) {
      plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
      text(.5,.5,paste("Heatmap cannot be computed",
                       "with less than two clusters.",sep="\n"))
    } else if (length(heatGenes()) < 1) {
      plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
      text(.5,.5,"There are no differentially expressed genes.")
    } else {
      if ("Unselected" %in% levels(clusts())) {
        tempLabRow <- c(paste(levels(clusts())[!levels(clusts()) == "Unselected"],
                              paste(sapply(switch(input$heatG,
                                                  deTissue=d$deTissue[[res()]],
                                                  deMarker=d$deMarker[[res()]],
                                                  deNeighb=deNeighb[[res()]]),nrow),"DE"),
                              sep=": "),"Unselected")
      } else {
        tempLabRow <- paste(paste0("Cluster ",levels(clusts())),
                            paste(sapply(switch(input$heatG,
                                                deTissue=d$deTissue[[res()]],
                                                deMarker=d$deMarker[[res()]],
                                                deNeighb=deNeighb[[res()]]),nrow),"DE"),
                            sep=": ")
      }
      heatmap.2(clustMeans(),Rowv=as.dendrogram(hC()),Colv=as.dendrogram(hG()),scale="column",
                margins=c(9,12),lhei=c(2,10),lwid=c(1,11),trace="none",
                keysize=1.5,density.info="none",key.par=list(mar=c(3,.5,2,.5),mgp=2:0),
                cexCol=1 + 1/log2(nrow(clustMeans())),cexRow=1 + 1/log2(ncol(clustMeans())),
                RowSideColors=clustCols(),labRow=tempLabRow,rowsep=sepClust(),col=viridis(100,d=-1))
    }
  }
  
  output$heatmap <- renderPlot({
    if (length(res()) > 0) {
      print(plot_heatmap())
    }
  })
  
  output$heatmapSave <- downloadHandler(
    filename="heatmap.pdf",
    content=function(file) {
      pdf(file,width=12,height=7)
      print(plot_heatmap())
      dev.off()
    }
  )
  
  output$deGeneSave <- downloadHandler(
    filename=function() { paste0(input$heatG,"_",input$DEclustNum,".txt") },
    content=function(file) {
      outTable <- switch(input$heatG,
                         deTissue=d$deTissue[[res()]][[input$DEclustNum]],
                         deMarker=d$deMarker[[res()]][[input$DEclustNum]],
                         deNeighb=deNeighb[[res()]][[input$DEclustNum]])
      write.table(outTable,file,quote=F,sep="\t",row.names=T,col.names=NA)
    }
  )
  

  #### clusterGenes ####
  output$genePlotClustSelect <- renderUI({
    if (length(res()) > 0) {
      selectInput("genePlotClust","Cluster:",selected=cSelected(),
                  choices=c("",levels(clusts())[!levels(clusts()) == "Unselected"]))
    }
  })

  cellMarkCols <- reactive(rainbow2(length(cellMarkers)))
  
  GOI <- eventReactive(input$GOIgo,grep(input$GOI,rownames(nge),value=T,ignore.case=T),ignoreNULL=F)
  
  plot_clusterGenes <- function() {
    doubleDot <- function(col1,col2) {
      upper.half.circle <- function(col1){  
        rs <- seq(0,pi,len=100) + pi/2
        xc <- 0+cos(rs) 
        yc <- 0+sin(rs) 
        polygon(xc,yc,col=col1,border=NA)
      } 
      lower.half.circle <- function(col2){ 
        rs <- seq(0,pi,len=100) + pi/2
        xc <- 0-cos(rs) 
        yc <- 0-sin(rs) 
        polygon(xc,yc,col=col2,border=NA)
      }
      upper.half.circle(col1)
      lower.half.circle(col2)
      rs <- seq(0,2*pi,len=200)
      polygon(cos(rs),sin(rs),border="white")
    }
    singleDot <- function(col1){  
      rs <- seq(0,2*pi,len=200)
      xc <- 0+cos(rs) 
      yc <- 0+sin(rs) 
      polygon(xc,yc,col=col1,border="white")
    } 
    par(mar=c(3,3,3,20),mgp=2:0)
    if (hiC() == "") {
      plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
      text(.5,.5,paste("Click a cell from a cluster on the tSNE plot above",
                       "to see gene expression for that cluster.",sep="\n"))
    } else {
      temp_ylab <- switch(as.character(exponent == exp(1)),
                          "TRUE"="(natural log scale)",
                          "FALSE"=paste0("(log",exponent," scale)"))
      plot(MDTC~DR,
           data=d$CGS[[res()]][[hiC()]][
             !((d$CGS[[res()]][[hiC()]]$cMu | d$CGS[[res()]][[hiC()]]$cMs) & 
                 d$CGS[[res()]][[hiC()]]$overCut),],
           col=alpha("black",0.3),
           xlab="Proportion of cells detecting gene",
           ylab=paste("Mean normalized gene expression of detected genes",temp_ylab))
      title(paste0("Cluster ", hiC(),": ",d$clusterID[[res()]][hiC()]),cex=1.2)
      mtext(paste("Cells:",sum(clusts()==hiC()),
                  "   Genes detected:",length(d$CGS[[res()]][[hiC()]]$DR)),side=3,line=0,cex=0.9)
      box(col=clustCols()[hiC()],lwd=2)

      if (input$cgLegend == "markers") {
        for (x in which(d$CGS[[res()]][[hiC()]]$cMu)) {
          my.symbols(x=d$CGS[[res()]][[hiC()]]$DR[x],
                     y=d$CGS[[res()]][[hiC()]]$MDTC[x],
                     symb=singleDot,inches=0.1,
                     MoreArgs=list(col1=cellMarkCols()[which(sapply(cellMarkersU,function(X) 
                       d$CGS[[res()]][[hiC()]]$genes[x] %in% X))]))
        }
        for (x in which(d$CGS[[res()]][[hiC()]]$cMs)) {
          temp <- unlist(strsplit(names(which(sapply(cellMarkersS,function(X) 
            d$CGS[[res()]][[hiC()]]$genes[x] %in% X))),"&"))
          my.symbols(x=d$CGS[[res()]][[hiC()]]$DR[x],
                     y=d$CGS[[res()]][[hiC()]]$MDTC[x],
                     symb=doubleDot,inches=0.1,
                     MoreArgs=list(col1=cellMarkCols()[as.integer(temp[1])],
                                   col2=cellMarkCols()[as.integer(temp[2])]))
        }
        for (x in which(d$CGS[[res()]][[hiC()]]$cMu & d$CGS[[res()]][[hiC()]]$overCut)) {
          text(x=d$CGS[[res()]][[hiC()]]$DR[x],y=d$CGS[[res()]][[hiC()]]$MDTC[x],
               labels=d$CGS[[res()]][[hiC()]]$genes[x],srt=315,cex=1.5,font=2,adj=c(1.1,-.1),
               col=cellMarkCols()[which(sapply(cellMarkersU,function(X) 
                 d$CGS[[res()]][[hiC()]]$genes[x] %in% X))])
        }
        for (x in which(d$CGS[[res()]][[hiC()]]$cMs & d$CGS[[res()]][[hiC()]]$overCut)) {
          text(x=d$CGS[[res()]][[hiC()]]$DR[x],y=d$CGS[[res()]][[hiC()]]$MDTC[x],
               labels=d$CGS[[res()]][[hiC()]]$genes[x],srt=315,cex=1.5,font=2,adj=c(1.1,-.1),
               col=cellMarkCols()[as.integer(temp[2])])
        }
        legend(x=1.05,y=max(d$CGS[[res()]][[hiC()]]$MDTC),xpd=NA,bty="n",ncol=1,
               pch=19,col=cellMarkCols(),legend=names(cellMarkersU))
        
      } else if (input$cgLegend == "heatmap") {
        degl <- rownames(d$CGS[[res()]][[hiC()]]) %in% 
          rownames(switch(input$heatG,
                          deTissue=d$deTissue[[res()]],
                          deMarker=d$deMarker[[res()]],
                          deNeighb=deNeighb[[res()]])[[hiC()]])[1:input$DEgeneCount]
        if (any(degl)) {
          points(x=d$CGS[[res()]][[hiC()]]$DR[degl],y=d$CGS[[res()]][[hiC()]]$MDTC[degl],
                 pch=16,cex=1.2,col="darkred")
          text(x=d$CGS[[res()]][[hiC()]]$DR[degl],y=d$CGS[[res()]][[hiC()]]$MDTC[degl],
               srt=315,cex=1.5,font=2,adj=c(1.1,-.1),col="darkred",
               labels=d$CGS[[res()]][[hiC()]]$genes[degl])
        }
        
      } else if (input$cgLegend == "regex" & length(GOI()) > 0) {
        degl <- which(rownames(nge) %in% GOI())
        points(x=d$CGS[[res()]][[hiC()]]$DR[degl],y=d$CGS[[res()]][[hiC()]]$MDTC[degl],
               pch=16,cex=1.2,col="darkred")
        text(x=d$CGS[[res()]][[hiC()]]$DR[degl],y=d$CGS[[res()]][[hiC()]]$MDTC[degl],
             srt=315,cex=1.5,font=2,adj=c(1.1,-.1),col="darkred",
             labels=d$CGS[[res()]][[hiC()]]$genes[degl])
      }
    }
  }
  
  output$clusterGenes <- renderPlot({
    if (length(res()) > 0) {
      print(plot_clusterGenes())
    }
  })
  
  output$clusterGenesSave <- downloadHandler(
    filename="clusterGenes.pdf",
    content=function(file) {
      pdf(file,width=12,height=7)
      print(plot_clusterGenes())
      dev.off()
    }
  )
  
  #### Gene Stats Plot ####
  cgGeneOpts <- reactive({
    t <- nearPoints(d$CGS[[res()]][[hiC()]],input$cgClick,xvar="DR",yvar="MDTC")
    return(t$genes)
  })
  
  output$cgSelect <- renderUI({
    if (length(res()) > 0) {
      if (input$boxplotGene == "click") {
        selectInput("cgGene",label="Gene:",choices=sort(cgGeneOpts()))
      } else if (input$boxplotGene == "regex") {
        selectInput("cgGene",label="Gene:",choices=sort(GOI()))
      }
    }
  })

  plot_geneTest <- function() {
    if (input$cgGene == "") {
      plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
      text(.5,.5,paste("Select a gene by either clicking on the plot above",
                       "or entering regular expression capturing your gene symbol of interest",
                       "then pick the gene from the list just above this figure",
                       "to see a comparison of that gene's expression across all clusters.",sep="\n"))
    } else {
      temp_ylab <- switch(as.character(exponent == exp(1)),
                          "TRUE"="(natural log scale)",
                          "FALSE"=paste0("(log",exponent," scale)"))
      temp_pos <- switch(as.character(length(levels(clusts())) > 1),"TRUE"=hC()$order,"FALSE"=1)
      layout(matrix(2:1,nrow=2),heights=c(1,4))
      par(mar=c(3,3,0,3),mgp=2:0)
      suppressWarnings(boxplot(vector("list",length(levels(clusts()))),
                               ylim=range(nge[input$cgGene,]),
                               ylab=paste(input$cgGene,"gene expression",temp_ylab),
                               xlab=NA,xaxt="n"))
      mtext(levels(clusts())[temp_pos],side=1,line=0,at=seq_along(temp_pos))
      mtext("Clusters, ordered by heatmap dendrogram",side=1,line=1)
      try(tempGeneName <- select(get(egDB),keys=input$cgGene,
                                 keytype="SYMBOL",column="GENENAME")$GENENAME,silent=T)
      if (exists("tempGeneName")) { 
        mtext(paste(paste("Gene name:",tempGeneName),collapse="\n"),
              side=1,line=2,font=2) 
      }
      if ("sct" %in% input$bxpOpts) {
        bxpCol <- alpha(clustCols(),.2)
      } else {
        bxpCol <- alpha(clustCols(),.8)
      }
      for (i in temp_pos) {
        boxplot(nge[input$cgGene,clusts() == levels(clusts())[i]],
                col=bxpCol[i],at=which(temp_pos == i),add=T,notch=T,outline=F)
        if ("sct" %in% input$bxpOpts) {
          points(jitter(rep(which(temp_pos == i),sum(clusts() == levels(clusts())[i])),amount=.2),
                 nge[input$cgGene,clusts() == levels(clusts())[i]],pch=20,col=alpha(clustCols()[i],.4))
        }
      }
      if ("rnk" %in% input$bxpOpts) {
        points(x=seq_along(d$CGS[[res()]]),
               y=sapply(d$CGS[[res()]][temp_pos],function(X) X[input$cgGene,"MTCrank"]) * 
                 max(nge[input$cgGene,]) + min(nge[input$cgGene,]),
               pch=25,cex=1.2,col="darkred",bg="firebrick2")
        axis(side=4,at=seq(0,1,.25) * max(nge[input$cgGene,]) + min(nge[input$cgGene,]),
             labels=percent(seq(0,1,.25)),col.ticks="darkred",col.axis="darkred")
        mtext(side=4,line=2,text="Quantile of gene expression per cluster",col="darkred")
      }
      if (length(temp_pos) > 1) { 
        par(new=F,mar=c(0,3,1,3))
        plot(as.dendrogram(hC()),leaflab="none") 
      }
    }
  }
  
  output$geneTest <- renderPlot({
    if (length(res()) > 0) {
      print(plot_geneTest())
    }
  })
  
  output$geneTestSave <- downloadHandler(
    filename="geneTest.pdf",
    content=function(file) {
      pdf(file,width=12,height=7)
      print(plot_geneTest())
      dev.off()
    }
  )
  
  
  ######## Distribution of genes of interest #########
  
  GOI1 <- eventReactive(input$GOI1go,
                        grep(input$GOI1,rownames(nge),value=T,ignore.case=T),
                        ignoreNULL=F)
  output$GOI1select <- renderUI({ 
    selectInput("goi1",label="Gene:",choices=sort(GOI1()),multiple=T)
  })
  
  GOI2 <- eventReactive(input$GOI2go,
                        grep(input$GOI2,rownames(nge),value=T,ignore.case=T),
                        ignoreNULL=F)
  output$GOI2select <- renderUI({ 
    selectInput("goi2",label="Gene:",choices=sort(GOI2()),multiple=T)
  })
  
  plot_goi <- function(goi) {
    if (length(goi) < 1) {
      plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
      text(.5,.5,paste("To search for your gene(s) of interest type a",
                       "search term (regex allowed) in the box above", 
                       "then select the gene(s) from the drop-down list",
                       "in the \"Gene:\" box above right.",sep="\n"))
    } else {
      if (length(goi) > 5) { goiL <- 5 } else { goiL <- length(goi) }
      if (goiL > 1) {
        gv <- apply(nge[goi,],2,max)
      } else {
        gv <- nge[goi,]
      }
      cv <- cut(gv,breaks=100,labels=F)
      par(mar=c(3,3,goiL+1,1),mgp=2:0)
      plot(dr_viz,pch=21,cex=1.3,xlab="tSNE_1",ylab="tSNE_2",
           col=viridis(100,.7,d=-1)[cv],bg=viridis(100,.3,d=-1)[cv])
      temp_yrange <- max(dr_viz[,2]) - min(dr_viz[,2])
      segments(x0=seq(quantile(range(dr_viz[,1]),.55),
                      quantile(range(dr_viz[,1]),.95),length.out=1000),
               y0=max(dr_viz[,2]) + temp_yrange * .045,
               y1=max(dr_viz[,2]) + temp_yrange * .065,
               col=viridis(1000,d=-1),xpd=NA)
      text(x=c(quantile(range(dr_viz[,1]),.55),
               quantile(range(dr_viz[,1]),.75),
               quantile(range(dr_viz[,1]),.95)),
           y=rep(max(dr_viz[,2]) + temp_yrange * .06,3),
           labels=c(round(min(gv),2),"Max expression per cell",round(max(gv),2)),pos=2:4,xpd=NA)
      try(tempGeneName <- 
            select(get(egDB),keys=goi,keytype="SYMBOL",column="GENENAME")$GENENAME,silent=T)
      if (exists("tempGeneName")) { 
        if (length(tempGeneName) > 4) { 
          tempGeneName[5] <- "and more..."; tempGeneName <- tempGeneName[1:5] 
        }
        title(paste(tempGeneName,collapse="\n"),line=0.25,adj=.01,font.main=1)
      }
    }
  }
  
  output$goiPlot1 <- renderPlot({
    if (input$plotClust1 == "clust" & length(res()) > 0) {
      print(plot_tsne())
      if (input$plotLabel1) { print(plot_tsne_labels()) }
    } else if (input$plotClust1 == "goi") {
      print(plot_goi(input$goi1))
      if (input$plotLabel1 & length(res()) > 0 & length(input$goi1) > 0) {
        print(plot_tsne_labels())
      }
    }
  })
  
  output$goiPlot1Save <- downloadHandler(
    filename="goi1.pdf",
    content=function(file) {
      pdf(file,width=7,height=7)
      if (input$plotClust1 == "clust" & length(res()) > 0) {
        print(plot_tsne())
        if (input$plotLabel1) { print(plot_tsne_labels()) }
      } else if (input$plotClust1 == "goi") {
        print(plot_goi(input$goi1))
        if (input$plotLabel1 & length(res()) > 0 & length(input$goi1) > 0) {
          print(plot_tsne_labels())
        }
      }
      dev.off()
    }
  )
  
  output$goiPlot2 <- renderPlot({
    if (input$plotClust2 == "clust" & length(res()) > 0) {
      print(plot_tsne())
      if (input$plotLabel2) { print(plot_tsne_labels()) }
    } else if (input$plotClust2 == "goi") {
      print(plot_goi(input$goi2))
      if (input$plotLabel2 & length(res()) > 0 & length(input$goi2) > 0) {
        print(plot_tsne_labels())
      }
    }
  })
  
  output$goiPlot2Save <- downloadHandler(
    filename="goi2.pdf",
    content=function(file) {
      pdf(file,width=7,height=7)
      if (input$plotClust2 == "clust" & length(res()) > 0) {
        print(plot_tsne())
        if (input$plotLabel2) { print(plot_tsne_labels()) }
      } else if (input$plotClust2 == "goi") {
        print(plot_goi(input$goi2))
        if (input$plotLabel2 & length(res()) > 0 & length(input$goi2) > 0) {
          print(plot_tsne_labels())
        }
      }
      dev.off()
    }
  )
  
  
  ######## Custom sets for DE #########
  selectedSets <- reactiveValues(a=NULL,b=NULL)
  
  plot_tsne_selDE <- function() {
    if (input$tsneSelDEcol == "") {
      id <- rep(1,nrow(md))
      idcol <- "grey20"
    } else if (is.factor(md[,input$tsneSelDEcol]) | is.character(md[,input$tsneSelDEcol])) {
      id <- as.factor(md[,input$tsneSelDEcol])
      if (length(levels(md[,input$tsneSelDEcol])) <= 8) {
        idcol <- brewer.pal(length(levels(md[,input$tsneSelDEcol])),
                            "Dark2")[1:length(levels(md[,input$tsneSelDEcol]))]
      } else {
        idcol <- rainbow2(length(levels(md[,input$tsneSelDEcol])))
      }
    } else {
      id <- cut(md[,input$tsneSelDEcol],100)
      idcol <- viridis(100,d=-1)
    }
    par(mar=c(3,3,3,1),mgp=2:0)
    plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
         xlim=range(dr_viz[,1]),ylim=range(dr_viz[,2]))
    points(dr_viz,pch=21,
           col=alpha(idcol,.8)[id],
           bg=alpha(idcol,0.4)[id])
    
    if (input$tsneSelDEcol == "") {
    } else if (is.factor(md[,input$tsneSelDEcol]) | is.character(md[,input$tsneSelDEcol])) {
      legend("topleft",bty="n",horiz=T,xpd=NA,inset=c(0,-.09),
             pch=21,col=idcol,pt.bg=alpha(idcol,0.5),
             title=input$tsneSelDEcol,legend=levels(md[,input$tsneSelDEcol]))
    } else {
      legend("topleft",bty="n",horiz=T,xpd=NA,inset=c(0,-.09),
             pch=21,col=viridis(3,d=-1),pt.bg=viridis(3,.5,d=-1),
             title=input$tsneSelDEcol,
             legend=c(round(min(md[,input$tsneSelDEcol]),2),
                      round((max(md[,input$tsneSelDEcol]) - 
                               min(md[,input$tsneSelDEcol])) / 2,2),
                      round(max(md[,input$tsneSelDEcol]),2)))
    }
    
    
    points(dr_viz[selectedSets$a,],pch=19,col="#a50026")
    points(dr_viz[selectedSets$b,],pch=19,col="#313695")
    points(dr_viz[intersect(selectedSets$a,selectedSets$b),],pch=19,col="#ffffbf")
    points(dr_viz[intersect(selectedSets$a,selectedSets$b),],pch=4,col="red")
    
    legend("topright",horiz=T,bty="n",xpd=NA,inset=c(0,-.09),
           title="Selected Cells",legend=c("Set A","Set B","Both"),
           pch=c(19,19,4),col=c("#a50026","#313695","red"))
  }
  output$tsneSelDE <- renderPlot({ print(plot_tsne_selDE()) })
  
  
  currSel <- reactive(rownames(brushedPoints(as.data.frame(dr_viz),
                                             input$tsneBrush,xvar="tSNE_1",yvar="tSNE_2")))
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
  output$textSetA <- renderText(paste(length(selectedSets$a),"cells in Set A."))
  output$textSetB <- renderText(paste(length(selectedSets$b),"cells in Set B."))
  output$textOverlap <- renderText(paste(length(intersect(selectedSets$a,selectedSets$b)),
                                         "cells in both sets.",
                                         "Cells must be assigned to a single set prior to calculation."))
  
  observeEvent(input$calcDE,{
    newRes <- paste0("Comp.",gsub("[^A-Za-z0-9]","",input$DEsetName))
    if (length(intersect(selectedSets$a,selectedSets$b)) > 0) {
      output$calcText <- renderText("Sets can't overlap (please assign red cells to only one set).")
    } else if (any(sapply(list(selectedSets$a,selectedSets$b),length) < 3)) {
      output$calcText <- renderText("Each set must contain at least 3 cells.")
    } else if (nchar(input$DEsetName) < 1) {
      output$calcText <- renderText("Please name this comparison (in text box above).")
    } else if (newRes %in% colnames(d$cl)) {
      output$calcText <- renderText("This comparison name has already been used.")
    } else {
      output$calcText <- renderText("")
      withProgress({
        temp <- rep("Unselected",nrow(d$cl))
        names(temp) <- rownames(d$cl)
        temp[selectedSets$a] <- "Set A"
        temp[selectedSets$b] <- "Set B"
        d$cl[[newRes]] <- factor(temp)
        
        #### Gene stats per set ####
        incProgress(amount=1/6,detail="Gene detection rate per set")
        setCells <- d$cl[,newRes] != "Unselected"
        DR <- apply(nge[,setCells],1,function(X) 
          tapply(X,d$cl[,newRes][setCells],function(Y) sum(Y>0)/length(Y)))
        
        incProgress(amount=1/6,detail="Mean detected gene expression per set")
        MDTC <- apply(nge[,setCells],1,function(X) 
          tapply(X,d$cl[,newRes][setCells],function(Y) {
            temp <- mean.logX(Y[Y>0])
            if (is.na(temp)) { temp <- 0 }
            return(temp)
          }))
        
        incProgress(amount=1/6,detail="Mean gene expression per set")
        MTC <- apply(nge,1,function(X) 
          tapply(X,d$cl[,newRes],mean.logX))
        
        d$CGS[[newRes]] <- sapply(levels(d$cl[,newRes])[1:2],function(X) 
          data.frame(DR=DR[X,],MDTC=MDTC[X,],MTC=MTC[X,]),simplify=F)
        d$CGS[[newRes]][["Unselected"]] <- data.frame(MTC=MTC["Unselected",])
        for (i in names(d$CGS[[newRes]])) {
          d$CGS[[newRes]][[i]]$MTCrank <- rank(d$CGS[[newRes]][[i]]$MTC,
                                               ties.method="min")/nrow(d$CGS[[newRes]][[i]])
          if (i == "Unselected") { next }
          d$CGS[[newRes]][[i]]$cMu <- rownames(d$CGS[[newRes]][[i]]) %in% unlist(cellMarkersU)
          d$CGS[[newRes]][[i]]$cMs <- rownames(d$CGS[[newRes]][[i]]) %in% unlist(cellMarkersS)
          d$CGS[[newRes]][[i]]$overCut <- d$CGS[[newRes]][[i]]$MTC > mean(d$CGS[[newRes]][[i]]$MTC)
          d$CGS[[newRes]][[i]]$genes <- rownames(d$CGS[[newRes]][[i]])
        }
        if (length(cellMarkers) < 1) {
          d$clusterID[[newRes]] <- rep("",nrow(cl))
        } else {
          d$clusterID[[newRes]] <- c(names(cellMarkers)[sapply(d$CGS[[newRes]][1:2],function(Y) 
            which.max(sapply(cellMarkers,function(X) median(Y$MTC[rownames(Y) %in% X]))))],
            "Unselected")
          names(d$clusterID[[newRes]]) <- names(d$CGS[[newRes]])
        }
        
        #### deTissue - DE per cluster vs all other data ####
        incProgress(amount=1/6,detail="DE vs tissue logGER calculations")
        deT_logGER <- sapply(levels(d$cl[,newRes])[1:2],function(i) 
          MTC[i,] - apply(nge[,d$cl[,newRes] != i],1,mean.logX))
        deT_genesUsed <- apply(deT_logGER,2,function(X) which(X > logGERthresh))  
        if (any(sapply(deT_genesUsed,length) < 1)) {
          stop(paste0("logGERthresh should be set to less than ",
                      min(apply(deT_logGER,2,function(X) max(abs(X)))),
                      ", the largest magnitude logGER between cluster ",
                      names(which.min(apply(deT_logGER,2,function(X) max(abs(X))))),
                      " and the remaining data."))
        }
        incProgress(amount=1/6,detail="DE vs tissue Wilcoxon rank sum calculations")
        deT_pVal <- sapply(levels(d$cl[,newRes])[1:2],function(i)
          apply(nge[deT_genesUsed[[i]],],1,function(X) 
            wilcox.test(X[d$cl[,newRes] == i],X[d$cl[,newRes] != i])$p.value),simplify=F)
        d$deTissue[[newRes]] <- sapply(levels(d$cl[,newRes])[1:2],function(i) 
          data.frame(logGER=deT_logGER[deT_genesUsed[[i]],i],
                     pVal=deT_pVal[[i]])[order(deT_pVal[[i]]),],simplify=F)
        tempQval <- tapply(
          p.adjust(do.call(rbind,d$deTissue[[newRes]])$pVal,"fdr"),
          rep(names(sapply(d$deTissue[[newRes]],nrow)),sapply(d$deTissue[[newRes]],nrow)),
          c)
        for (i in names(d$deTissue[[newRes]])) { 
          d$deTissue[[newRes]][[i]] <- d$deTissue[[newRes]][[i]][tempQval[[i]] <= WRSTalpha,]
          d$deTissue[[newRes]][[i]]$qVal <- tempQval[[i]][tempQval[[i]] <= WRSTalpha] 
        }
        
        #### deMarker - DE per cluster vs each other cluster #### 
        incProgress(amount=1/6,detail="Calculating Set A vs Set B")

        deM_dDR <- DR["Set A",] - DR["Set B",]
        deM_logGER <- MTC["Set A",] - MTC["Set B",]
        deM_genesUsed <- switch(threshType,
                                dDR=which(abs(deM_dDR) > dDRthresh),
                                logGER=which(abs(deM_logGER) > logGERthresh))
        if (length(deM_genesUsed) < 1) {
          stop("Gene filtering threshold is set too high.")
        }

        deM_pVal <- apply(nge[deM_genesUsed,],1,function(X) 
          wilcox.test(X[d$cl[,newRes] == "Set A"],
                      X[d$cl[,newRes] == "Set B"])$p.value)
        
        temp_deVS <- data.frame(dDR=deM_dDR[deM_genesUsed],
                                logGER=deM_logGER[deM_genesUsed],
                                pVal=deM_pVal)[order(deM_pVal),]
        temp_deVS$qVal <- p.adjust(temp_deVS$pVal,"fdr")

        d$deMarker[[newRes]] <- list(
          "Set A"=temp_deVS[temp_deVS[,threshType] > 0 & temp_deVS$qVal <= WRSTalpha,],
          "Set B"=temp_deVS[temp_deVS[,threshType] < 0 & temp_deVS$qVal <= WRSTalpha,]
        )
        d$deMarker[[newRes]][["Set B"]]$dDR <- d$deMarker[[newRes]][["Set B"]]$dDR * -1
        d$deMarker[[newRes]][["Set B"]]$logGER <- d$deMarker[[newRes]][["Set B"]]$logGER * -1

        selectedSets$a <- selectedSets$b <- NULL
      },message="DE calculations:")      

      
    }
  })
  observeEvent(input$updateForViz, {
    withProgress({
      new_cl <- d$cl[input$res]
      new_CGS <- list()
      for (i in names(d$CGS[[input$res]])) {
        new_CGS[[input$res]][[i]] <- 
          d$CGS[[input$res]][[i]][colnames(d$CGS[[input$res]][[i]]) %in% c("DR","MDTC","MTC")]
      }
      new_deTissue <- d$deTissue[input$res]
      new_deMarker <- d$deMarker[input$res]
      incProgress(.5)
      save(new_cl,new_CGS,new_deTissue,new_deMarker,
           file=paste0(dataPath,dataTitle,"_selDE_",sub("Comp.","",input$res,fixed=T),".RData"))
    },message=paste0(
      "Saving ",dataTitle,"_selDE_",sub("Comp.","",input$res,fixed=T),".RData to ",dataPath))
  })
  
}


########## ShinyApp ##########
shinyApp(ui = ui, server = server)