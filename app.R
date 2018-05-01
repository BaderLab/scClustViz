
#### UI ####
ui <- fixedPage(
  fixedRow(
    titlePanel(dataTitle)),
  hr(),
  
  fixedRow(titlePanel("Cluster Resolution Selection"),
    column(6,radioButtons("res","Resolution:",choices=names(deMarker),inline=T,selected=savedRes),
           fixedRow(column(6,actionButton("go","View clusters at this resolution"),align="left"),
                    column(6,actionButton("save","Save this resolution as default"),align="right")),
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
    
  fixedRow(
    titlePanel("Cell-type Clusters"),
    column(6,downloadButton("libSizeSave","Save as PDF"),
           plotOutput("libSize",height="550px"),align="right")#,
#    column(6,downloadButton("libSizeSave","Save as PDF"),
#           plotOutput("libSize",height="550px"),align="right")
  ),
  
  fixedRow(
    column(6,fixedRow(
      if (length(cellMarkers) > 0) {
        column(4,radioButtons("tsneLabels","Labels:",inline=T,
                              choices=list("Cluster numbers"="cn","Cluster annotations"="ca")))
      } else {
        column(4,radioButtons("tsneLabels","Labels:",inline=T,
                              choices=list("Cluster numbers"="cn")))
      },
      column(2,downloadButton("tsneSave","Save as PDF"),align="right")
    ),
    plotOutput("tsne",height="600px",click="tsneClick")),
    column(6,align="right",
           fixedRow(downloadButton("cellCycleSave","Save as PDF"),
                    plotOutput("cellCycle",height="600px"),align="right"))
  ),
  hr(),

  fixedRow(
    titlePanel("Cluster-wise Gene Stats"),
    column(5,fixedRow(
      if (length(cellMarkers) > 0) {
        column(2,radioButtons("cgLegend",label="Highlighted genes:",
                              choices=c("Cell-type markers"="markers",
                                        "Top DE genes (from heatmap)"="heatmap",
                                        "Gene symbol (regex)"="regex")))
      } else {
        column(2,radioButtons("cgLegend",label="Highlighted genes:",
                              choices=c("Top DE genes (from heatmap)"="heatmap",
                                        "Gene symbol (regex)"="regex")))
      },
      column(3,align="right",
             textInput("GOI","Gene symbol (regex)"),
             downloadButton("clusterGenesSave","Save as PDF"))
    ),
    plotOutput("clusterGenes",height="500px",click="cgClick")
    ),
    column(7,fixedRow(
      column(2,radioButtons("boxplotGene",label="Gene of interest:",
                            choices=c("Click from plot on left"="click",
                                      "Gene symbol (regex)"="regex")),
             downloadButton("geneTestSave","Save as PDF")),
      column(5,uiOutput("cgRadio"))
    ),
    plotOutput("geneTest",height="500px"))
  ),
  hr(),
  
  fixedRow(
    column(2,radioButtons("heatG","Heapmap Genes:",
                          choices=list("DE vs tissue average"="deTissue",
                                       "Marker genes"="deMarker"))),
    column(2,uiOutput("DEclustSelect"),align="right"),
    column(2,downloadButton("deGeneSave","Download gene list"),
           downloadButton("heatmapSave","Save as PDF"),align="right"), 
    column(6,uiOutput("DEgeneSlider"))
  ),
  fixedRow(plotOutput("heatmap",height="600px")),
  hr(),
  
  fixedRow(
    titlePanel("Genes of Interest"),
    fixedRow(6,fixedRow(
      column(2,align="left",
             checkboxInput("plotClust1",label="Plot gene expression overlay",value=F),
             checkboxInput("plotLabel1",label="Include labels from tSNE plot above",value=T),
             downloadButton("goiPlot1Save","Save as PDF")),
      column(4,textInput("GOI1",label="Gene symbol (regex):"),align="right")
    ),
    plotOutput("goiPlot1",height="700px")
    ),
    column(6,fixedRow(
      column(2,align="left",
             checkboxInput("plotClust2",label="Plot gene expression overlay",value=F),
             checkboxInput("plotLabel2",label="Include labels from tSNE plot above",value=T),
             downloadButton("goiPlot2Save","Save as PDF")),
      column(4,textInput("GOI2",label="Gene symbol (regex):"),align="right")
    ),
    plotOutput("goiPlot2",height="700px")
    )
  )
)


#### Server ####
server <- function(input,output,session) {
  output$dimension_display <- renderText({
    paste(input$dimension[1], input$dimension[2], input$dimension[2]/input$dimension[1])
  })
  
  clustCols <- reactive({      
    if (length(levels(cl[,input$res])) <= 8) {
      brewer.pal(length(levels(cl[,input$res])),"Dark2")[1:length(levels(cl[,input$res]))]
    } else {
      rainbow2(length(levels(cl[,input$res])))
    }
  })
  
  #### cqPlot ####
  plot_cqPlot <- function() {
    numClust <- sapply(cl,function(X) length(levels(X)))
    numDEgenes <- lapply(get(input$deType),function(X) sapply(X,nrow))
    x <- which(names(numClust) == input$res)
    toplim <- c(21,max(unlist(numDEgenes)) + 20)
    botlim <- c(-1,21)
    
    layout(matrix(c(3,3,1,2),2),widths=c(1,29))
    par(mar=c(0.2,2,1,1),mgp=2:0)
    plot(x=numClust,y=sapply(numDEgenes,median),type="l",
         xlim=range(numClust)+c(-.5,.5),ylim=toplim,yaxs="i",xaxt="n",ylab=NA)
    abline(h=seq(0,max(unlist(numDEgenes)),10),lty=3,col=alpha(1,0.3))
    for (i in names(numDEgenes)[-x]) {
      boxplot(numDEgenes[[i]],add=T,at=numClust[[i]])
    }
    boxplot(numDEgenes[[x]],add=T,at=numClust[[x]],
            border="red")
    
    par(mar=c(3,2,0.2,1),mgp=2:0)
    plot(x=numClust,y=sapply(numDEgenes,median),type="l",
         xlim=range(numClust)+c(-.5,.5),ylim=botlim,yaxs="i",xlab="Number of clusters",ylab=NA)
    abline(h=seq(0,max(unlist(numDEgenes)),10),lty=3,col=alpha(1,0.3))
    for (i in names(numDEgenes)[-x]) {
      boxplot(numDEgenes[[i]],add=T,at=numClust[[i]])
    }
    boxplot(numDEgenes[[x]],add=T,at=numClust[[x]],
            border="red")
    
    par(mar=c(3,0,1,0))
    plot.new()
    mtext(switch(input$deType,
                 "deMarker"="Positive DE genes per cluster to all other clusters",
                 "deNeighb"="Positive DE genes per cluster to nearest cluster")
          ,side=2,line=-1.5)
  }
  
  output$cqPlot <- renderPlot({
    print(plot_cqPlot())
  })
  
  output$cqPlotSave <- downloadHandler(
    filename="cqPlot.pdf",
    content=function(file) {
      pdf(file,width=12,height=9)
      print(plot_cqPlot())
      dev.off()
    }
  )
  
  #### sil ####
  plot_sil <- function() {
    tempSil <- silhouette(as.integer(cl[,input$res]),dist=silDist)
    par(mar=c(4,0,2,1),mgp=2:0)
    plot(tempSil,beside=T,border=NA,main=NA,col=clustCols(),do.n.k=T)
  }
  
  output$sil <- renderPlot({
    print(plot_sil())
  })
  
  output$silSave <- downloadHandler(
    filename="sil.pdf",
    content=function(file) {
      pdf(file,width=9,height=12)
      print(plot_sil())
      dev.off()
    }
  )
  
  #### res buttons ####
  res <- eventReactive(input$go,input$res)
  
  observeEvent(input$save,{
    savedRes <<- input$res #<<- updates variable outside scope of function (ie. global environment)
    save(savedRes,file=paste0(dataPath,dataTitle,"_savedRes.RData"))
  })
  
  #### Cluster Data ####
  clusts <- reactive(cl[,res()])
  
  #### tSNE ####
  transp <- reactiveValues(bg=rep(0.5,ncol(nge)),col=rep(1,ncol(nge)))
  
  plot_tsne_labels <- function() {
    if (input$tsneLabels == "ca") {
      temp <- sapply(unique(clusterID[[res()]]),function(X) names(which(clusterID[[res()]] == X)))
      text(apply(dr_viz,2,function(Y) tapply(Y,apply(sapply(temp,function(X) clusts() %in% X),1,which),mean)),
        labels=names(temp),font=2,cex=1.5)
    } else if (input$tsneLabels == "cn") {
      text(apply(dr_viz,2,function(X) tapply(X,clusts(),mean)),
           labels=levels(clusts()),font=2,cex=1.5)
    } else {
      legend("center",legend="You changed the label choice names...")
    }
  }
  
  plot_tsne <- function() {
    par(mar=c(4,3,3,1),mgp=2:0)
    plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
         main=paste("tSNE at",res(),"using",ncol(dr_clust),"PCs"),
         xlim=range(dr_viz[,1]),ylim=range(dr_viz[,2]))
    points(dr_viz,pch=21,
           col=alpha(clustCols()[clusts()],transp$col),
           bg=alpha(clustCols()[clusts()],transp$bg))
    plot_tsne_labels()
    if (!is.na(hiC())) {
      mtext(side=3,line=-1,text=paste("Cluster",hiC(),"-",
                                      clusterID[[res()]][hiC()],"-",
                                      sum(clusts() == hiC()),"cells"))
    }
  }
  
  output$tsne <- renderPlot({
    print(plot_tsne())
  })
  
  output$tsneSave <- downloadHandler(
    filename="tsne.pdf",
    content=function(file) {
      pdf(file,width=10,height=10)
      print(plot_tsne())
      dev.off()
    }
  )
  
  #### tsneClick ####
  tsneClick <- reactiveValues(cl=NULL)
  
  observeEvent(input$tsneClick,{ tsneClick$cl <- input$tsneClick })
  
  hiC <- reactive({
    t <- nearPoints(as.data.frame(dr_viz),tsneClick$cl,xvar="tSNE_1",yvar="tSNE_2",threshold=5)
    t2 <- cl[rownames(t)[1],res()]
    return(t2)
  })
  
  ci <- reactive({
    if (is.na(hiC())) {
      rep(F,length(clusts()))
    } else {
      clusts() == hiC()
    }
  })
  
  observeEvent(input$tsneClick,{
    transp$bg <- 0.5
    transp$cl <- 1
    if (any(ci())) {
      transp$bg[!ci()] <- 0.1
      transp$col[!ci()] <- 0.2
    }
  })
  
  #### Cell cycle tSNE ####
  plot_cellCycle <- function() {
    layout(cbind(2:1),heights=c(1,7))
    par(mar=c(3,3,0,1),mgp=2:0)
    plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
         xlim=range(dr_viz[,1]),ylim=range(dr_viz[,2]))
    if (any(ci())) {
      points(dr_viz[!ci(),],pch=21,
             col=viridis(3,.1)[c(3,1,2)][cycScores$phases[!ci()]],
             bg=viridis(3,0.05)[c(3,1,2)][cycScores$phases[!ci()]])
      points(dr_viz[ci(),],pch=21,
             col=viridis(3,.8)[c(3,1,2)][cycScores$phases[ci()]],
             bg=viridis(3,0.4)[c(3,1,2)][cycScores$phases[ci()]])
    } else {
      points(dr_viz,pch=21,
             col=viridis(3,.8)[c(3,1,2)][cycScores$phases],
             bg=viridis(3,0.4)[c(3,1,2)][cycScores$phases])
    }
    par(mar=c(0,0,0,0))
    plot.new()
    legend("top",bty="n",horiz=T,pch=c(NA,21,21,21),
           legend=c("Cell cycle\nphase:",levels(cycScores$phases)),
           col=c(NA,viridis(3)[c(3,1,2)]),pt.bg=c(NA,viridis(3,0.5)[c(3,1,2)]))
  }
  
  output$cellCycle <- renderPlot({
    print(plot_cellCycle())
  })
  
  output$cellCycleSave <- downloadHandler(
    filename="cellCycle.pdf",
    content=function(file) {
      pdf(file,width=10,height=10)
      print(plot_cellCycle())
      dev.off()
    }
  )
  
  #### Library Size Plot ####
  densGene <- reactive({
    if (any(ci())) {
      list(density(md[!ci(),"total_features"]),
           density(md[ci(),"total_features"]))
    } else {
      list(density(md$total_features))
    }
  })
  
  densUMI <- reactive({
    if (any(ci())) {
      list(density(md[!ci(),"total_counts"]),
           density(md[ci(),"total_counts"]))
    } else {
      list(density(md$total_counts))
    }
  })
  
  plot_libSize <- function() {
    layout(matrix(c(2,1,0,3),2),c(6,2),c(2,6))
    par(mar=c(3,3,0,0),mgp=2:0)
    plot(total_features~total_counts,data=md[!ci(),],
         pch=21,col=alpha("black",0.2),bg=alpha("black",0.1),cex=1.2,
         xlab="Library Size",ylab="Genes Detected")
    points(total_features~total_counts,data=md[ci(),],
           pch=21,col=alpha("red",0.4),bg=alpha("red",0.2),cex=1.2)
    if (any(ci())) {
      legend("topleft",bty="n",pch=21,col="red",pt.bg=alpha("red",0.5),
             legend=paste("Cluster",hiC(),"-",clusterID[[res()]][hiC()]))
    }
    par(mar=c(0,3,1,0))
    plot(x=NULL,y=NULL,ylab="Density",xaxt="n",
         xlim=range(md$total_counts),
         ylim=c(min(sapply(densUMI(),function(X) min(X$y))),
                max(sapply(densUMI(),function(X) max(X$y)))))
    for (x in 1:length(densUMI())) {
      lines(densUMI()[[x]],col=c("black","red")[x],lwd=3)
    }
    par(mar=c(3,0,0,1))
    plot(x=NULL,y=NULL,xlab="Density",yaxt="n",
         xlim=c(min(sapply(densGene(),function(X) min(X$y))),
                max(sapply(densGene(),function(X) max(X$y)))),
         ylim=range(md$total_features))
    for (x in 1:length(densGene())) {
      lines(x=densGene()[[x]]$y,y=densGene()[[x]]$x,col=c("black","red")[x],lwd=3)
    }
  }
  
  output$libSize <- renderPlot({
    print(plot_libSize())
  })
  
  output$libSizeSave <- downloadHandler(
    filename="libSize.pdf",
    content=function(file) {
      pdf(file,width=10,height=10)
      print(plot_libSize())
      dev.off()
    }
  )
  
  #### clusterGenes ####
  cellMarkCols <- reactive(gg_colour_hue(length(cellMarkers)))
  
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
    with(CGS()[[hiC()]],{
      plot(x=DR[!((cMu | cMs) & overCut)],
           y=MDTC[!((cMu | cMs) & overCut)],
           col=alpha("black",0.3),
           xlab="Proportion of cells detecting gene",
           ylab="Mean normalized gene expression of detected genes")
      title(paste0("Cluster ", hiC(),": ",clusterID[[res()]][hiC()]),cex=1.2)
      mtext(paste("Cells:",sum(clusts()==hiC()),
                  "   Genes detected:",length(CGS()[[hiC()]]$DR)),side=3,line=0,cex=0.9)
      box(col=clustCols()[hiC()],lwd=2)
      lines(x=mean(MTC)/seq(min(MDTC),max(MDTC)*.75,by=.01),
            y=seq(min(MDTC),max(MDTC)*.75,by=.01),
            lty=2,lwd=2,col=alpha("red",0.5))
      
      if (input$cgLegend == "markers") {
        for (x in which(cMu)) {
          my.symbols(x=DR[x],y=MDTC[x],symb=singleDot,inches=0.1,
                     MoreArgs=list(col1=cellMarkCols()[which(sapply(cellMarkersU,function(X) genes[x] %in% X))]))
        }
        for (x in which(cMs)) {
          temp <- unlist(strsplit(names(which(sapply(cellMarkersS,function(X) genes[x] %in% X))),"&"))
          my.symbols(x=DR[x],y=MDTC[x],symb=doubleDot,inches=0.1,
                     MoreArgs=list(col1=cellMarkCols()[as.integer(temp[1])],
                                   col2=cellMarkCols()[as.integer(temp[2])]))
        }
        for (x in which(cMu & overCut)) {
          text(x=DR[x],y=MDTC[x],labels=genes[x],
               srt=315,cex=1.5,font=2,adj=c(1.1,-.1),
               col=cellMarkCols()[which(sapply(cellMarkersU,function(X) genes[x] %in% X))])
        }
        for (x in which(cMs & overCut)) {
          text(x=DR[x],y=MDTC[x],labels=genes[x],
               srt=315,cex=1.5,font=2,adj=c(1.1,-.1),
               col=cellMarkCols()[as.integer(temp[2])])
        }
        legend("top",inset=.05,bty="n",ncol=2,#floor((length(cellMarkCols())+3)/4),
               pch=c(rep(19,times=length(cellMarkCols())),NA,NA),
               lty=c(rep(NA,times=length(cellMarkCols())),2,NA),
               lwd=c(rep(NA,times=length(cellMarkCols())),2,NA),
               col=c(cellMarkCols(),"red",NA),
               legend=c(names(cellMarkersU),"Mean of mean expression",
                        paste(sum(CGS()[[hiC()]]$overCut),"over threshold")))
      } 
      
      else if (input$cgLegend == "heatmap") {
        degl <- rownames(CGS()[[hiC()]]) %in% 
          switch(input$heatG,
                 deTissue=if (ncol(deTissue()[[hiC()]]) == 0) { NA } else { 
                   deTissue()[[hiC()]][1:input$DEgeneCount,"gene"] },
                 deMarker=if (ncol(deMarker()[[hiC()]]) == 0) { NA } else { 
                   deMarker()[[hiC()]][1:input$DEgeneCount,"gene"] })
        if (any(degl)) {
          points(x=DR[degl],y=MDTC[degl],
                 pch=16,cex=1.2,col="darkred")
          text(x=DR[degl],y=MDTC[degl],
               srt=315,cex=1.5,font=2,adj=c(1.1,-.1),
               col="darkred",labels=genes[degl])
        }
        legend("top",inset=.05,bty="n",horiz=T,
               lty=c(2,NA),lwd=c(2,NA),col=c("red",NA),
               legend=c("Mean of mean expression",
                        paste(sum(CGS()[[hiC()]]$overCut),"over threshold")))
      } 
      
      else if (input$cgLegend == "regex"){
        degl <- grep(input$GOI,genes)
        points(x=DR[degl],y=MDTC[degl],
               pch=16,cex=1.2,col="darkred")
        text(x=DR[degl],y=MDTC[degl],
             srt=315,cex=1.5,font=2,adj=c(1.1,-.1),
             col="darkred",labels=genes[degl])
        legend("top",inset=.05,bty="n",horiz=T,
               lty=c(2,NA),lwd=c(2,NA),col=c("red",NA),
               legend=c("Mean of mean expression",
                        paste(sum(CGS()[[hiC()]]$overCut),"over threshold")))
      } 
      
      else {
        legend("center",legend="You changed the choice names...")
      }
    })
  }
  
  output$clusterGenes <- renderPlot({
    print(plot_clusterGenes())
  })
  
  output$clusterGenesSave <- downloadHandler(
    filename="clusterGenes.pdf",
    content=function(file) {
      pdf(file,width=12,height=9)
      print(plot_clusterGenes())
      dev.off()
    }
  )
  
  #### Gene Stats Plot ####
  cgGeneOpts <- reactive({
    t <- nearPoints(CGS()[[hiC()]],input$cgClick,xvar="DR",yvar="MDTC")
    return(t$genes)
  })
  
  output$cgRadio <- renderUI({
    if (input$boxplotGene == "click") {
      radioButtons("cgGene",label="Gene:",choices=sort(cgGeneOpts()),inline=T)
    }
  })
  
  plot_geneTest <- function() {
    if (input$boxplotGene == "click") {
      goi <- input$cgGene
    } else if (input$boxplotGene == "regex") {
      goi <- grep(input$GOI,rownames(nge),value=T)
    }
    if (length(goi) > 5) { goiL <- 5 } else { goiL <- length(goi) }
    
    temp <- sapply(CGS(), function(X) {
      if (goi %in% X$genes) {
        which(X[order(X$MTC,decreasing=T),]$genes == goi)
      } else {
        nrow(X)
      }
    })
    names(temp) <- levels(clusts())
    
    layout(matrix(2:1,nrow=2),heights=c(1,4))
    par(mar=c(2+goiL,3,0,3),mgp=2:0)
    if (goiL > 1) {
      plot(x=NULL,y=NULL,xlim=c(1,length(temp)),
           ylim=range(nge[goi,]),
           ylab="Max gene expression per cell",xlab=NA,xaxt="n")
    } else {
      plot(x=NULL,y=NULL,xlim=c(1,length(temp)),
           ylim=range(nge[goi,]),
           ylab=paste(goi,"gene expression (log2)"),xlab=NA,xaxt="n")
    }
    mtext(levels(clusts())[hC()$order],side=1,line=0,at=seq_along(temp))
    mtext("Clusters, ordered by heatmap dendrogram",side=1,line=1)
    try(tempGeneName <- select(get(egDB),keys=goi,keytype="SYMBOL",column="GENENAME")$GENENAME,silent=T)
    if (exists("tempGeneName")) { 
      if (length(tempGeneName) > 5) { tempGeneName[5] <- "and more..."; tempGeneName <- tempGeneName[1:5] }
      mtext(paste(paste("Gene name:",tempGeneName),collapse="\n"),
            side=1,line=goiL+1,font=2) 
    }
    for (i in hC()$order) {
      if (goiL > 1) {
        boxplot(apply(nge[goi,clusts() == levels(clusts())[i]],2,max),
                at=which(hC()$order == i),add=T,col=clustCols()[i])
      } else {
        boxplot(nge[goi,clusts() == levels(clusts())[i]],
                at=which(hC()$order == i),add=T,col=clustCols()[i])
      }
    }
    par(new=T)
    plot(x=seq_along(CGS()),y=sapply(CGS()[hC()$order],function(X) max(X[goi,"MTCrank"])),
         axes=F,xlab=NA,ylab=NA,ylim=0:1,pch=25,cex=1.2,col="darkred",bg="coral")
    axis(side=4,col.ticks="darkred",col.axis="darkred")
    mtext(side=4,line=2,text="Gene expression quantile per cluster",col="darkred")
    
    par(new=F,mar=c(0,3,1,3))
    plot(as.dendrogram(hC()),leaflab="none",xaxs="i")
  }
  
  output$geneTest <- renderPlot({
    print(plot_geneTest())
  })
  
  output$geneTestSave <- downloadHandler(
    filename="geneTest.pdf",
    content=function(file) {
      pdf(file,width=12,height=9)
      print(plot_geneTest())
      dev.off()
    }
  )
  
  
  #### Heatmap genes ####
  deTissue <- reactive({ get(load(paste0(dataPath,"precalc_",gsub(".","",res(),fixed=T),"_deTissue.RData"))) })
  #deMarker <- reactive({ get(load(paste0(dataPath,"precalc_",gsub(".","",res(),fixed=T),"_deMarker.RData"))) })
  
  output$DEgeneSlider <- renderUI({ 
    switch(input$heatG,
           deTissue=sliderInput("DEgeneCount",min=2,max=max(sapply(deTissue(),nrow)),value=5,step=1,ticks=T,width="100%",
                                label=HTML(paste("Positive differential gene expression of cluster over tissue",
                                                 "# of genes per cluster to show",sep="<br/>"))),
           deMarker=sliderInput("DEgeneCount",min=2,max=max(sapply(deMarker(),nrow)),value=5,step=1,ticks=T,width="100%",
                                label=HTML(paste("Positive differential gene expression between cluster and all other clusters",
                                                 "# of genes per cluster to show",sep="<br/>"))))
  })
  
  output$DEclustSelect <- renderUI({
    selectInput("DEclustNum","Cluster # for gene list",choices=seq_along(deTissue()),selectize=F)
  })
  
  heatGenes <- reactive({
    temp <- switch(input$heatG,
                   deTissue=unique(unlist(lapply(deTissue(),function(X) 
                     if (ncol(X) == 0) { NA } else { X[1:input$DEgeneCount,"gene"] }))),
                   deMarker=unique(unlist(lapply(deMarker(),function(X) 
                     if (ncol(X) == 0) { NA } else { X[1:input$DEgeneCount,"gene"] }))))
    temp <- temp[!is.na(temp)]
    return(temp)
  })
  
  clustMeans <- reactive({ #This only works if input is in ascending order of adjusted p value.
    temp <- sapply(CGS(),function(X) X[heatGenes(),"MTC"])
    rownames(temp) <- heatGenes()
    return(t(temp))
  })
  
  hC <- reactive(hclust(dist(clustMeans()),"single"))
  hG <- reactive(hclust(dist(t(clustMeans())),"complete"))
  
  sepClust <- reactive({
    if (is.na(hiC())) {
      return(c(NA,NA))
    } else {
      return(nrow(clustMeans()) - 
               c(which(hC()$order == hiC()) - 1,
                 which(hC()$order == hiC())))
    }
  })
  
  plot_heatmap <- function() {
    tempLabRow <- switch(input$heatG,
                         deTissue=paste(paste0("#",seq_along(deTissue())),
                                        paste(sapply(deTissue(),function(X) 
                                          if (ncol(X) > 0) { nrow(X) } else { return(0) }),"DE"),sep=": "),
                         deMarker=paste(paste0("#",seq_along(deMarker())),
                                        paste(sapply(deMarker(),function(X) 
                                          if (ncol(X) > 0) { nrow(X) } else { return(0) }),"DE"),sep=": "))
    heatmap.2(clustMeans(),Rowv=as.dendrogram(hC()),Colv=as.dendrogram(hG()),scale="column",
              col="viridis",trace="none",margins=c(9,7),keysize=1,lhei=c(2,10),lwid=c(1,11),
              cexCol=1 + 1/log2(nrow(clustMeans())),cexRow=1 + 1/log2(ncol(clustMeans())),
              RowSideColors=clustCols(),labRow=tempLabRow,rowsep=sepClust())
  }
  
  output$heatmap <- renderPlot({
    print(plot_heatmap())
  })
  
  output$heatmapSave <- downloadHandler(
    filename="heatmap.pdf",
    content=function(file) {
      pdf(file,width=9,height=12)
      print(plot_heatmap())
      dev.off()
    }
  )
  
  output$deGeneSave <- downloadHandler(
    filename=function() { paste0(input$heatG,"_",input$DEclustNum,".txt") },
    content=function(file) {
      outTable <- switch(input$heatG,
                         deTissue=deTissue()[[as.integer(input$DEclustNum)]],
                         deMarker=deMarker()[[as.integer(input$DEclustNum)]])
      write.table(outTable,file,quote=F,sep="\t",row.names=T,col.names=NA)
    }
  )
  
  
  #### Cluster Explorer ####
  goi1 <- reactive(input$GOI1)
  goi2 <- reactive(input$GOI2)
  
  plot_goi1 <- function() {
    if (!input$plotClust1) {
      print(plot_tsne())
    } else {
      goi <- grep(goi1(),rownames(nge),value=T)
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
      if (input$plotLabel1) { plot_tsne_labels() }
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
      try(tempGeneName <- select(get(egDB),keys=goi,keytype="SYMBOL",column="GENENAME")$GENENAME,silent=T)
      if (exists("tempGeneName")) { 
        if (length(tempGeneName) > 5) { tempGeneName[5] <- "and more..."; tempGeneName <- tempGeneName[1:5] }
        title(paste(tempGeneName,collapse="\n"),line=0.25,adj=.01,font.main=1)
      }
    }
  }
  
  output$goiPlot1 <- renderPlot({
    print(plot_goi1())
  })
  
  output$goiPlot1Save <- downloadHandler(
    filename="goi1.pdf",
    content=function(file) {
      pdf(file,width=10,height=10)
      print(plot_goi1())
      dev.off()
    }
  )
  
  plot_goi2 <- function() {
    if (!input$plotClust2) {
      print(plot_tsne())
    } else {
      goi <- grep(goi2(),rownames(nge),value=T)
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
      if (input$plotLabel2) { plot_tsne_labels() }
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
      try(tempGeneName <- select(get(egDB),keys=goi,keytype="SYMBOL",column="GENENAME")$GENENAME,silent=T)
      if (exists("tempGeneName")) { 
        if (length(tempGeneName) > 5) { tempGeneName[5] <- "and more..."; tempGeneName <- tempGeneName[1:5] }
        title(paste(tempGeneName,collapse="\n"),line=0.25,adj=.01,font.main=1)
      }
    }
  }
  
  output$goiPlot2 <- renderPlot({
    print(plot_goi2())
  })
  
  output$goiPlot2Save <- downloadHandler(
    filename="goi2.pdf",
    content=function(file) {
      pdf(file,width=10,height=10)
      print(plot_goi2())
      dev.off()
    }
  )
  
}

shinyApp(ui=ui,server=server)


#### TESTING ####
# Put this in the UI:
#  fluidRow(textOutput("testing")),#### TESTING ####
#  hr(),#### TESTING
# Put this in the Server at the point you want to check on.
#output$testing <- renderPrint(input$tsneClick)  #### TESTING ####

