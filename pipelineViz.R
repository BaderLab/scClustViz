dataName <- "yourData"
setwd(paste0("../",dataName,"/output")) # point to output directory

### Set species of data
#species <- "human"
species <- "mouse"

library(shiny)
library(cluster)
library(gplots)
library(viridis)
library(RColorBrewer)
library(scales)
library(TeachingDemos)
library(vioplot)
if (species == "human") {
  library(org.Hs.eg.db) # from Bioconductor (for human data)
  egDB <- "org.Hs.eg.db"
} else if (species == "mouse") {
  library(org.Mm.eg.db) # from Bioconductor (for mouse data)
  egDB <- "org.Mm.eg.db"
} else { } #Gene name lookup won't work.
library(Seurat) #see http://satijalab.org/seurat/install.html

if (!exists("eb1S")) {
  load(paste0(dataName,"_eb1S.RData"))
}

if (!exists("cycScores")) {
  load(paste0(dataName,"_cycScores.RData"))
}

### Generate list of cell-type markers. 
## List should be character vectors of gene symbols 
## (check to see that you're using the right symbol - HGNC/MGI symbols),
## with list names representing the cell type.  The list below is an example:
cellMarkers <- list("Cortical precursors"=c("Mki67","Sox2","Pax6","Pcna","Nes","Cux1","Cux2"),
                    "Interneurons"=c("Gad1","Gad2","Npy","Sst","Lhx6","Tubb3","Rbfox3","Dcx"),
                    "Cajal-Retzius neurons"="Reln",
                    "Intermediate progenitors"="Eomes",
                    "Projection neurons"=c("Tbr1","Satb2","Fezf2","Bcl11b","Tle4",
                                           "Nes","Cux1","Cux2","Tubb3","Rbfox3","Dcx"),
                    "Oligodendrocyte precursors"=c("Cspg4","Olig2","Pdgfra"),
                    "Oligodendrocytes"=c("Mbp","Mog","Plp1","Mag"),
                    "Astrocytes"=c("Aldh1l1","Gfap","Slc1a3","Glul"),
                    "Microglia"="Cx3cr1")

### Handling of non-unique markers - visualizations exist for genes that mark up to two different cell types.
## Genes marking three or more cell types may break this or the gene visualization
## (although I may be able to bump it to 4 cell types if necessary).
cellMarkersS <- apply(combn(seq_along(cellMarkers),2),2,function(X) do.call(intersect,unname(cellMarkers[X])))
names(cellMarkersS) <- apply(combn(seq_along(cellMarkers),2),2,function(X) paste(X,collapse="&"))
cellMarkersS <- cellMarkersS[sapply(cellMarkersS,length) > 0]
cellMarkersU <- lapply(cellMarkers,function(X) X[!X %in% unlist(cellMarkersS)])

gg_colour_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

if (file.exists(paste0(dataName,"_savedRes.RData"))) {
  load(paste0(dataName,"_savedRes.RData"))
} else {
  savedRes <- NULL
}

shinyApp(
  #### UI ####
  ui = fluidPage(
    fluidRow(
      titlePanel(dataName) # Could also be a string with a title of your choice.
    ),
    hr(),
    fluidRow(
      column(6,
             titlePanel("Cluster Resolution Selection"),
             radioButtons("res","Resolution:",choices=names(minDEgenes),inline=T,selected=savedRes),
             fluidRow(column(6,actionButton("go","View clusters at this resolution")),
                      column(6,actionButton("save","Save this resolution as default"))
             ),
             plotOutput("cqPlot",height="400px")),
      column(6,plotOutput("sil",height="600px"))
    ),
    hr(),
    fluidRow(
      titlePanel("Cell-type Clusters"),
      column(6,
             radioButtons("tsneLabels","Labels:",inline=T,
                          choices=list("Cluster numbers"="cn","Cluster annotations"="ca")),
             plotOutput("tsne",height="700px",click="tsneClick"),
             fluidRow(
               column(6,plotOutput("cellCycle",height="350px")),
               column(6,plotOutput("libSize",height="350px"))
             )),
      column(6,
             fluidRow(
               column(6,radioButtons("heatG","Heapmap Genes:",
                                     choices=list("DE vs tissue"="deGall",
                                                  "DE vs all other cell types"="deGvs"),
                                     inline=T)),
               column(6,uiOutput("DEgeneSlider"))
             ),
             plotOutput("heatmap",height="1000px"))
    ),
    hr(),
    fluidRow(
      titlePanel("Cluster-wise Gene Stats"),
      column(5,
             fluidRow(
               column(6,radioButtons("cgLegend",label="Highlighted genes:",
                                     choices=c("Cell-type markers",
                                               "Top DE genes (from heatmap)",
                                               "GOI (symbols,)"))),
               column(6,textInput("GOI","GOI (symbols,)"))
             ),
             plotOutput("clusterGenes",height="500px",click="cgClick")
      ),
      column(7,
             uiOutput("cgRadio"),
             plotOutput("geneTest",height="500px"))
    ),
    hr(),
    fluidRow(
      titlePanel("Genes of Interest"),
      column(6,
             fluidRow(
               column(4,checkboxInput("plotClust1",label="Clusters?",value=T)),
               column(8,textInput("GOI1",label="Gene of interest (1):"))
             ),
             plotOutput("goiPlot1",height="700px")
      ),
      column(6,
             fluidRow(
               column(4,checkboxInput("plotClust2",label="Clusters?",value=T)),
               column(8,textInput("GOI2",label="Gene of interest (1):"))
             ),
             plotOutput("goiPlot2",height="700px")
      )
    )
  ),
  #### Server ####
  server = function(input,output,session) {
    
    #### cqPlot ####
    numClust <- apply(eb1S@data.info[,grepl("res",colnames(eb1S@data.info))],2,function(Y) length(unique(Y)))
    numDEgenes <- lapply(uniqDE,function(X) sapply(X,length))
    output$cqPlot <- renderPlot({
      x <- which(names(numClust) == input$res)
      par(mar=c(3,3,1,1),mgp=2:0)
      plot(x=numClust,y=sapply(numDEgenes,median),type="l",
           xlim=range(numClust),ylim=range(unlist(numDEgenes)),
           xlab="Number of clusters",
           ylab="DE genes (@ 1% FDR) per cluster to all other clusters")
      abline(h=seq(0,max(unlist(numDEgenes)),10),lty=3,col=alpha(1,0.3))
      for (i in names(numDEgenes)[-x]) {
        boxplot(numDEgenes[[i]],add=T,at=numClust[[i]])
      }
      boxplot(numDEgenes[[x]],add=T,at=numClust[[x]],
              border="red")
    })
    
    #### res buttons ####
    res <- eventReactive(input$go,input$res)
    
    observeEvent(input$save,{
      savedRes <<- input$res
      save(savedRes,file=paste0(dataName,"_savedRes.RData"))
    })
    
    #### Cluster Data ####
    clusts <- reactive(eb1S@data.info[,res()])
    CGS <- reactive({
      load(paste0(dataName,"_precalc_",gsub(".","",res(),fixed=T),"_CGS.RData"))
      for (i in seq_along(CGS)) {
        CGS[[i]]$MTCrank <- rank(CGS[[i]]$MTC,ties.method="min")/nrow(CGS[[i]])
        CGS[[i]]$cMu <- rownames(CGS[[i]]) %in% unlist(cellMarkersU)
        CGS[[i]]$cMs <- rownames(CGS[[i]]) %in% unlist(cellMarkersS)
        CGS[[i]]$overCut <- CGS[[i]]$MTC > mean(CGS[[i]]$MTC)
        CGS[[i]]$genes <- rownames(CGS[[i]])
      }
      return(CGS)
    })
    output$test <- renderPrint({
      sapply(CGS(),dim)
    })
    
    clusterID <- reactive({
      names(sapply(CGS(),function(X) which.max(sapply(cellMarkers,function(Y) median(X$MTC[X$genes %in% Y])))))
    })
    clustCols <- reactive({      
      if (length(levels(clusts())) <= 8) {
        brewer.pal(length(levels(clusts())),"Dark2")
      } else {
        gg_colour_hue(length(levels(clusts())))
      }
    })
    
    #### sil ####
    output$sil <- renderPlot({
      tempDist <- dist(eb1S@pca.rot[,seq(1,maxPCt)],method="euclidean")
      tempSil <- silhouette(as.integer(clusts()),dist=tempDist)
      par(mar=c(4,0,2,1),mgp=2:0)
      plot(tempSil,beside=T,border=NA,main=NA,col=clustCols(),do.n.k=T)
    })
    
    #### tSNE ####
    transp <- reactiveValues(bg=rep(0.5,ncol(eb1S@data)),col=rep(1,ncol(eb1S@data)))
    output$tsne <- renderPlot({
      par(mar=c(4,3,3,1),mgp=2:0)
      plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
           main=paste(paste("tSNE of",dataName,"at",res()),
                      paste(maxPCt,"PCs from",dim(eb1S@pca.x)[1],"genes"),
                      sep="\n"),
           xlim=range(eb1S@tsne.rot[,1]),ylim=range(eb1S@tsne.rot[,2]))
      points(eb1S@tsne.rot,pch=21,
             col=alpha(clustCols()[clusts()],transp$col),
             bg=alpha(clustCols()[clusts()],transp$bg))
      if (input$tsneLabels == "ca") {
        temp <- lapply(unique(clusterID()),function(X) which(clusterID() == X))
        names(temp) <- unique(clusterID())
        text(x=tapply(FUN=mean,X=eb1S@tsne.rot[,1],
                      INDEX=apply(sapply(temp,function(X) clusts() %in% X),1,which)),
             y=tapply(FUN=mean,X=eb1S@tsne.rot[,2],
                      INDEX=apply(sapply(temp,function(X) clusts() %in% X),1,which)),
             labels=names(temp),col="black",font=2,cex=1.5)
      } else if (input$tsneLabels == "cn") {
        text(x=tapply(FUN=mean,X=eb1S@tsne.rot[,1],INDEX=clusts()),
             y=tapply(FUN=mean,X=eb1S@tsne.rot[,2],INDEX=clusts()),
             labels=levels(clusts()),col="black",font=2,cex=1.5)
      } else {
        legend("center",legend="You changed the choice names...")
      }
      if (!is.na(hiC())) {
        mtext(side=3,line=-1,text=paste("Cluster",hiC(),"-",
                                        clusterID()[hiC()],"-",
                                        sum(clusts() == hiC()),"cells"))
      }
    })
    
    #### tsneClick ####
    tsneClick <- reactiveValues(cl=NULL)
    observeEvent(input$tsneClick,{ tsneClick$cl <- input$tsneClick })
    hiC <- reactive({
      t <- nearPoints(eb1S@tsne.rot,tsneClick$cl,xvar="tSNE_1",yvar="tSNE_2",threshold=5)
      t2 <- eb1S@data.info[rownames(t)[1],res()]
      return(as.integer(as.character(t2)))
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
    
    #### Context Plots ####
    output$cellCycle <- renderPlot({
      layout(cbind(2:1),heights=c(1,7))
      par(mar=c(3,3,0,1),mgp=2:0)
      plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
           xlim=range(eb1S@tsne.rot[,1]),ylim=range(eb1S@tsne.rot[,2]))
      if (any(ci())) {
        points(eb1S@tsne.rot[!ci(),],pch=21,
               col=viridis(3,.1)[c(3,1,2)][cycScores$phases[!ci()]],
               bg=viridis(3,0.05)[c(3,1,2)][cycScores$phases[!ci()]])
        points(eb1S@tsne.rot[ci(),],pch=21,
               col=viridis(3,.8)[c(3,1,2)][cycScores$phases[ci()]],
               bg=viridis(3,0.4)[c(3,1,2)][cycScores$phases[ci()]])
      } else {
        points(eb1S@tsne.rot,pch=21,
               col=viridis(3,.8)[c(3,1,2)][cycScores$phases],
               bg=viridis(3,0.4)[c(3,1,2)][cycScores$phases])
      }
      par(mar=c(0,3,0,1))
      plot.new()
      legend("top",bty="n",horiz=T,pch=c(NA,21,21,21),
             legend=c("Cell cycle\nphase:",levels(cycScores$phases)),
             col=c(NA,viridis(3)[c(3,1,2)]),pt.bg=c(NA,viridis(3,0.5)[c(3,1,2)]))
    })
    
    densGene <- reactive({
      if (any(ci())) {
        list(density(pDat[!ci(),"total_features"]),
             density(pDat[ci(),"total_features"]))
      } else {
        list(density(pDat$total_features))
      }
    })
    densUMI <- reactive({
      if (any(ci())) {
        list(density(pDat[!ci(),"total_counts"]),
             density(pDat[ci(),"total_counts"]))
      } else {
        list(density(pDat$total_counts))
      }
    })
    output$libSize <- renderPlot({
      layout(matrix(c(2,1,0,3),2),c(6,2),c(2,6))
      par(mar=c(3,3,0,0),mgp=2:0)
      plot(total_features~total_counts,data=pDat[!ci(),],
           pch=21,col=alpha("black",0.2),bg=alpha("black",0.1),cex=1.2,
           xlab="Library Size",ylab="Genes Detected")
      points(total_features~total_counts,data=pDat[ci(),],
             pch=21,col=alpha("red",0.4),bg=alpha("red",0.2),cex=1.2)
      if (any(ci())) {
        legend("topleft",bty="n",pch=21,col="red",pt.bg=alpha("red",0.5),
               legend=paste("Cluster",hiC(),"-",clusterID()[hiC()]))
      }
      par(mar=c(0,3,1,0))
      plot(x=NULL,y=NULL,ylab="Density",xaxt="n",
           xlim=range(pDat$total_counts),
           ylim=c(min(sapply(densUMI(),function(X) min(X$y))),
                  max(sapply(densUMI(),function(X) max(X$y)))))
      for (x in 1:length(densUMI())) {
        lines(densUMI()[[x]],col=c("black","red")[x],lwd=3)
      }
      par(mar=c(3,0,0,1))
      plot(x=NULL,y=NULL,xlab="Density",yaxt="n",
           xlim=c(min(sapply(densGene(),function(X) min(X$y))),
                  max(sapply(densGene(),function(X) max(X$y)))),
           ylim=range(pDat$total_features))
      for (x in 1:length(densGene())) {
        lines(x=densGene()[[x]]$y,y=densGene()[[x]]$x,col=c("black","red")[x],lwd=3)
      }
    })
    
    #### Heatmap genes ####
    heatType <- reactive({ input$heatG == "deGall" })
    deG <- reactive({
      if (heatType()) {
        return(get(load(paste0(dataName,"_precalc_",gsub(".","",res(),fixed=T),"_deGall.RData"))))
      } else {
        return(get(load(paste0(dataName,"_precalc_",gsub(".","",res(),fixed=T),"_deGvs.RData"))))
      }
    })
    output$DEgeneSlider <- renderUI({
      if (heatType()) {
        sliderInput("DEgeneCount",label="Top differentially expressed genes per cluster",
                    min=2,max=max(sapply(deG(),nrow)),value=5,step=1,ticks=FALSE)
      } else {
        sliderInput("DEgeneCount",label="Top differentially expressed genes per cluster",
                    min=2,max=max(sapply(uniqDE[[res()]],length)),value=5,step=1,ticks=FALSE)
      }
    })
    heatGenes <- reactive({
      n <- input$DEgeneCount
      if (heatType()){
        return(lapply(deG(),function(X) rownames(X[order(X$fdr),])[1:n]))
      } else {
        temp <- lapply(names(uniqDE[[res()]])[order(as.integer(names(uniqDE[[res()]])))],
                       function(cl) sapply(uniqDE[[res()]][[cl]], 
                                           function(X) max(deG()[deG()$posClust == cl & deG()$gene == X,"fdr"])))
        temp2 <- lapply(temp,function(X) if(length(X)>0) {names(sort(X))[1:n]})
        return(lapply(temp2,function(X) X[!is.na(X)]))
      }
    })
    clustMeans <- reactive({
      temp <- sapply(CGS(),function(X) X$MTC[X$genes %in% unique(unlist(heatGenes()))])
      rownames(temp) <- CGS()[[1]]$genes[CGS()[[1]]$genes %in% unique(unlist(heatGenes()))]
      return(temp)
    })
    hG <- reactive(hclust(dist(clustMeans()),"complete"))
    hC <- reactive(hclust(dist(t(clustMeans())),"single"))
    
    sepClust <- reactive({
      if (is.na(hiC())) {
        return(c(NA,NA))
      } else {
        return(c(which(hC()$order == hiC()) - 1,
                 which(hC()$order == hiC())))
      }
    })
    
    output$heatmap <- renderPlot({
      if (heatType()) {
        tempLabCol <- paste(paste0("#",seq_along(deG())),
                            paste(sapply(deG(),nrow),"DE"),sep=": ")
      } else {
        tempDE <- uniqDE[[res()]][order(as.integer(names(uniqDE[[res()]])))]
        tempLabCol <- paste(paste0("#",names(tempDE)),
                            paste(sapply(tempDE,length),"DE"),sep=": ")
      }
      heatmap.2(clustMeans(),Rowv=as.dendrogram(hG()),Colv=as.dendrogram(hC()),scale="row",
                col="viridis",trace="none",margins=c(8,8),keysize=1,lhei=c(1,9),
                cexRow=1 + 1/log2(nrow(clustMeans())),cexCol=1 + 1/log2(ncol(clustMeans())),
                ColSideColors=clustCols(),labCol=tempLabCol,colsep=sepClust())
    })
    
    #### clusterGenes ####
    cellMarkCols <- reactive(gg_colour_hue(length(cellMarkers)))
    
    output$clusterGenes <- renderPlot({
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
        title(paste0("Cluster ", hiC(),": ",clusterID()[hiC()]),cex=1.2)
        mtext(paste("Cells:",sum(clusts()==hiC()),
                    "   Genes detected:",length(CGS()[[hiC()]]$DR)),side=3,line=0,cex=0.9)
        box(col=clustCols()[hiC()],lwd=2)
        lines(x=mean(MTC)/seq(min(MDTC),max(MDTC)*.75,by=.01),
              y=seq(min(MDTC),max(MDTC)*.75,by=.01),
              lty=2,lwd=2,col=alpha("red",0.5))
        if (input$cgLegend == "Cell-type markers") {
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
        } else if (input$cgLegend == "Top DE genes (from heatmap)") {
          degl <- rownames(CGS()[[hiC()]]) %in% heatGenes()[[hiC()]]
          points(x=DR[degl],y=MDTC[degl],
                 pch=16,cex=1.2,col="darkred")
          text(x=DR[degl],y=MDTC[degl],
               srt=315,cex=1.5,font=2,adj=c(1.1,-.1),
               col="darkred",labels=genes[degl])
          legend("top",inset=.05,bty="n",horiz=T,
                 lty=c(2,NA),lwd=c(2,NA),col=c("red",NA),
                 legend=c("Mean of mean expression",
                          paste(sum(CGS()[[hiC()]]$overCut),"over threshold")))
        } else if (input$cgLegend == "GOI (symbols,)"){
          degl <- sapply(strsplit(input$GOI,split=", ?",perl=T)[[1]],
                         function(X) which(X == genes))
          points(x=DR[degl],y=MDTC[degl],
                 pch=16,cex=1.2,col="darkred")
          text(x=DR[degl],y=MDTC[degl],
               srt=315,cex=1.5,font=2,adj=c(1.1,-.1),
               col="darkred",labels=genes[degl])
          legend("top",inset=.05,bty="n",horiz=T,
                 lty=c(2,NA),lwd=c(2,NA),col=c("red",NA),
                 legend=c("Mean of mean expression",
                          paste(sum(CGS()[[hiC()]]$overCut),"over threshold")))
        } else {
          legend("center",legend="You changed the choice names...")
        }
      })
    })
    
    #### Gene Stats Plot ####
    cgGeneOpts <- reactive({
      t <- nearPoints(CGS()[[hiC()]],input$cgClick,xvar="DR",yvar="MDTC")
      return(t$genes)
    })
    
    output$cgRadio <- renderUI({
      radioButtons("cgGene",label="Gene:",choices=sort(cgGeneOpts()),inline=T)
    })
    
    output$geneTest <- renderPlot({
      temp <- sapply(CGS(), function(X) {
        if (input$cgGene %in% X$genes) {
          which(X[order(X$MTC,decreasing=T),]$genes == input$cgGene)
        } else {
          nrow(X)
        }
      })
      names(temp) <- levels(clusts())
      
      layout(matrix(2:1,nrow=2),heights=c(1,3))
      par(mar=c(3,3,0,3),mgp=2:0)
      plot(x=NULL,y=NULL,xlim=c(1,length(temp)),
           ylim=range(eb1S@data[input$cgGene,]),
           ylab=paste(input$cgGene,"gene expression"),xlab=NA,xaxt="n")
      mtext(hC()$order,side=1,line=0,at=seq_along(temp))
      mtext("Clusters, ordered by heatmap dendrogram",side=1,line=1)
      try(tempGeneName <- select(get(egDB),keys=input$cgGene,keytype="SYMBOL",column="GENENAME")$GENENAME,silent=T)
      if (exists("tempGeneName")) { mtext(paste("Gene name:",tempGeneName),side=1,line=2,font=2) }
      for (i in hC()$order) {
        vioplot(eb1S@data[input$cgGene,clusts() == i],
                at=which(hC()$order == i),add=T,col=clustCols()[i])
      }
      par(new=T)
      plot(x=seq_along(CGS()),y=sapply(CGS()[hC()$order],function(X) X[input$cgGene,"MTCrank"]),
           axes=F,xlab=NA,ylab=NA,ylim=0:1,pch=25,cex=1.2,col="black",bg="darkred")
      axis(side=4,col.ticks="darkred",col.axis="darkred")
      mtext(side=4,line=2,text="Gene expression quantile per cluster",col="darkred")
      
      par(new=F,mar=c(0,3,1,3))
      plot(as.dendrogram(hC()),leaflab="none",xaxs="i")
    })
    
    #### Cluster Explorer ####
    goi1 <- reactive(input$GOI1)
    output$goiPlot1 <- renderPlot({
      if (input$plotClust1) {
        par(mar=c(4,3,3,1),mgp=2:0)
        plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
             main=paste(maxPCt,"PCs from",dim(eb1S@pca.x)[1],"genes"),
             xlim=range(eb1S@tsne.rot[,1]),ylim=range(eb1S@tsne.rot[,2]))
        points(eb1S@tsne.rot,pch=21,
               col=clustCols()[clusts()],
               bg=alpha(clustCols()[clusts()],0.5))
        
        temp <- lapply(unique(clusterID()),function(X) which(clusterID() == X))
        names(temp) <- unique(clusterID())
        text(x=tapply(FUN=mean,X=eb1S@tsne.rot[,1],
                      INDEX=apply(sapply(temp,function(X) clusts() %in% X),1,which)),
             y=tapply(FUN=mean,X=eb1S@tsne.rot[,2],
                      INDEX=apply(sapply(temp,function(X) clusts() %in% X),1,which)),
             labels=names(temp),col="black",font=2,cex=1.5)
      } else {
        iH <- 101 - cut(eb1S@data[goi1(),],breaks = 100,labels=F)
        plot(eb1S@tsne.rot[order(iH,decreasing=T),],xlab=NA,ylab=NA,pch=21,
             col=viridis(100,0.5)[sort(iH,decreasing=T)],
             bg=viridis(100,0.3)[sort(iH,decreasing=T)])
      }
    })
    
    goi2 <- reactive(input$GOI2)
    output$goiPlot2 <- renderPlot({
      if (input$plotClust2) {
        par(mar=c(4,3,3,1),mgp=2:0)
        plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
             main=paste(maxPCt,"PCs from",dim(eb1S@pca.x)[1],"genes"),
             xlim=range(eb1S@tsne.rot[,1]),ylim=range(eb1S@tsne.rot[,2]))
        points(eb1S@tsne.rot,pch=21,
               col=clustCols()[clusts()],
               bg=alpha(clustCols()[clusts()],0.5))
        
        temp <- lapply(unique(clusterID()),function(X) which(clusterID() == X))
        names(temp) <- unique(clusterID())
        text(x=tapply(FUN=mean,X=eb1S@tsne.rot[,1],
                      INDEX=apply(sapply(temp,function(X) clusts() %in% X),1,which)),
             y=tapply(FUN=mean,X=eb1S@tsne.rot[,2],
                      INDEX=apply(sapply(temp,function(X) clusts() %in% X),1,which)),
             labels=names(temp),col="black",font=2,cex=1.5)
      } else {
        iH <- 101 - cut(eb1S@data[goi2(),],breaks = 100,labels=F)
        plot(eb1S@tsne.rot[order(iH,decreasing=T),],xlab=NA,ylab=NA,pch=21,
             col=viridis(100,0.5)[sort(iH,decreasing=T)],
             bg=viridis(100,0.3)[sort(iH,decreasing=T)])
      }
    })
    
  }
)
