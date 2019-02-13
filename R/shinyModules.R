#' @include deTest.R speadLabels.R
NULL


# Cluster Solution DE boxplots -------------
plot_clustSep <- function(sCVdL,DEtype,FDRthresh,res,Xlim,Ylim) {
  numClust <- sapply(sCVdL,function(X) length(levels(Clusters(X))))
  if (is.null(Xlim)) { Xlim <- range(numClust) }
  bpData <- sapply(sCVdL,function(X) switch(DEtype,
                                            # DR=DEdist(X,"DR"),
                                            # MGE=DEdist(X,"MGE"),
                                            # PCA=DEdist(X,getEmb(inD,Param(sCVdL[[1]],"DRforClust"))),
                                            # scoreDE=as.vector(as.dist(DEdist(X))),
                                            DEneighb=sapply(DEneighb(X,FDRthresh),nrow),
                                            DEmarker=sapply(DEmarker(X,FDRthresh),nrow),
                                            silWidth=Silhouette(X)[,"sil_width"]),
                   simplify=F)
  if (is.null(Ylim)) { Ylim <- range(unlist(bpData)) }
  
  if (grepl("^Comp",res)) {
    par(mar=c(3,3,2,1))
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,paste("Press 'View clusters at this resolution'",
                     "to view the comparison",
                     sub("Comp.","",res,fixed=T),sep="\n"))
  } else {
    par(mar=c(3,3,2,1),mgp=2:0)
    if (DEtype == "silWidth") {
      plot(x=NA,y=NA,xlim=Xlim + c(-.5,.5),ylim=Ylim,
           xlab="Number of clusters",ylab="Silhouette width per cluster")
    } else {
      plot(x=numClust,y=sapply(bpData,median),type="l",
           xlim=Xlim + c(-.5,.5),ylim=Ylim,xlab="Number of clusters",
           ylab=switch(DEtype,
                       # DR="Distance between clusters by gene detection rates",
                       # MGE="Distance between clusters by mean gene expression",
                       scoreDE="Distance between clusters by differential expression score",
                       DEmarker="Positive DE genes per cluster to all other clusters",
                       DEneighb="Positive DE genes per cluster to nearest cluster"))
    }
    abline(h=seq(0,max(unlist(bpData)),switch(as.character(diff(Ylim) > 1000),
                                              "FALSE"=10,"TRUE"=100)),
           lty=3,col=alpha(1,0.3))
    for (i in names(bpData)[names(bpData) != res]) {
      boxplot(bpData[[i]],add=T,at=numClust[i],yaxt="n")
    }
    if (any(names(bpData) == res)) {
      if (DEtype == "silWidth") {
        boxplot(bpData[[res]],add=T,at=numClust[res],border="red")
      } else {
        boxplot(bpData[[res]],add=T,at=numClust[res],border="red",outline=F)
        points(jitter(rep(numClust[res],length(bpData[[res]])),amount=.2),
               bpData[[res]],col=alpha("red",.5),pch=20)
      }
    }
    if (DEtype == "silWidth") {
      temp_avSil <- sapply(bpData,mean)
      lines(numClust,y=temp_avSil,type="b",col="darkred",pch=16)
      points(numClust[res],temp_avSil[res],col="red",pch=16)
      legend(x=par("usr")[2],y=par("usr")[4],
             xjust=1,yjust=0.2,xpd=NA,bty="n",horiz=T,
             legend=c("Average silhouette width",
                      paste("Selected resolution:",res)),
             col=c("darkred","red"),pch=c(16,0),lty=c(1,NA))
    } else {
      legend(x=par("usr")[2],y=par("usr")[4],
             xjust=1,yjust=0.2,xpd=NA,bty="n",
             legend=paste("Selected resolution:",res),
             col="red",pch=0)
    }
  }
}

# plot_clustSep(sCVdL,DEtype="deMarker",FDRthresh=0.05,res="res.0.8",Ylim=NULL,Xlim=NULL)

# Silhouette plot ------
plot_sil <- function(sCVd) {
  par(mar=c(4.5,.5,1.5,1.5),mgp=2:0)
  plot(Silhouette(sCVd),
       beside=T,border=NA,main=NA,
       col=rainbow2(length(levels(Clusters(sCVd)))),
       do.n.k=T)
}


# tsnePlot -------------------
tsne_labels <- function(sCVd,cell_coord,lab_type) {
  if (lab_type == "ClusterNames") {
    temp_labelNames <- sapply(unique(attr(Clusters(sCVd),"ClusterNames")),function(X) 
      names(which(attr(Clusters(sCVd),"ClusterNames") == X)),simplify=F)
    temp_labels <- apply(cell_coord,2,function(Y) 
      tapply(Y,apply(sapply(temp_labelNames,function(X) Clusters(sCVd) %in% X),1,which),mean))
    if (!is.matrix(temp_labels)) { temp_labels <- rbind(temp_labels) }
    rownames(temp_labels) <- names(temp_labelNames)
  } else if (lab_type == "ClusterNamesAll") {
    temp_labels <- apply(cell_coord,2,function(X) tapply(X,Clusters(sCVd),mean))
    if (!is.matrix(temp_labels)) { temp_labels <- rbind(temp_labels) }
    rownames(temp_labels) <- attr(Clusters(sCVd),"ClusterNames")
  } else if (lab_type == "Clusters") {
    temp_labels <- apply(cell_coord,2,function(X) tapply(X,Clusters(sCVd),mean))
    if (!is.matrix(temp_labels)) { temp_labels <- rbind(temp_labels) }
    rownames(temp_labels) <- levels(Clusters(sCVd))
  } else {
    stop("lab_type should be one of 'ClusterNames', 'ClusterNamesAll', or 'Clusters'.")
  }
  return(temp_labels)
}

plot_tsne <- function(cell_coord,md,md_title,md_log,label,
                      sel_cells,sel_cells_A,sel_cells_B) {
  if (is.null(md_title)) {
    id <- as.factor(md)
    idcol <- rainbow2(length(levels(id)))
    if (any(is.na(id))) {
      levels(id) <- c(levels(id),"Unselected")
      id[is.na(id)] <- "Unselected"
      idcol <- c(idcol,"grey80")
    }
    par(mar=c(3,3,2,1),mgp=2:0)
  } else if (is.factor(md) | is.character(md)) {
    id <- as.factor(md)
    if (length(levels(md)) <= 8) {
      idcol <- RColorBrewer::brewer.pal(length(levels(id)),"Dark2")[1:length(levels(id))]
    } else {
      idcol <- rainbow2(length(levels(id)))
    }
    par(mar=c(3,3,ceiling(length(levels(id))/4)+.5,1),mgp=2:0)
  } else {
    if (md_log) {
      id <- cut(log10(md),100)
    } else {
      id <- cut(md,100)
    }
    idcol <- viridis::viridis(100,d=-1)
    par(mar=c(3,3,2,1),mgp=2:0)
  }
  if (missing(sel_cells)) { sel_cells <- character() }
  
  plot(x=NULL,y=NULL,xlab="tSNE_1",ylab="tSNE_2",
       xlim=range(cell_coord[,1]),ylim=range(cell_coord[,2]))
  if (length(sel_cells) > 0) {
    points(cell_coord[!rownames(cell_coord) %in% sel_cells,],pch=21,
           col=alpha(idcol,.6)[id[!rownames(cell_coord) %in% sel_cells]],
           bg=alpha(idcol,0.3)[id[!rownames(cell_coord) %in% sel_cells]])
    points(cell_coord[sel_cells,],pch=21,cex=1.3,
           col=alpha(idcol,1)[id[rownames(cell_coord) %in% sel_cells]],
           bg=alpha(idcol,0.6)[id[rownames(cell_coord) %in% sel_cells]])
  } else {
    points(cell_coord,pch=21,col=alpha(idcol,.8)[id],bg=alpha(idcol,0.4)[id])
  }
  
  if (!missing(sel_cells_A) & !missing(sel_cells_B)) {
    points(x=cell_coord[sel_cells_A,1],
           y=cell_coord[sel_cells_A,2],
           pch=19,col="#a50026")
    points(x=cell_coord[sel_cells_B,1],
           y=cell_coord[sel_cells_B,2],
           pch=19,col="#313695")
    points(x=cell_coord[intersect(sel_cells_A,sel_cells_B),1],
           y=cell_coord[intersect(sel_cells_A,sel_cells_B),2],
           pch=19,col="#ffffbf")
    points(x=cell_coord[intersect(sel_cells_A,sel_cells_B),1],
           y=cell_coord[intersect(sel_cells_A,sel_cells_B),2],
           pch=4,col="red")
  }
  if (!is.null(label)) {
    text(label,labels=rownames(label),font=2,cex=1.2)
  }
  if (is.null(md_title)) {
  } else if (is.factor(md) | is.character(md)) {
    legend(x=par("usr")[2],y=par("usr")[4],
           xjust=1,yjust=0.2,xpd=NA,bty="n",
           ncol=switch(as.character(length(levels(md)) < 4),"TRUE"=length(levels(md)),"FALSE"=4),
           legend=levels(id),pch=21,col=idcol,pt.bg=alpha(idcol,0.5))
    mtext(md_title,side=3,adj=0,font=2,line=ceiling(length(levels(id))/4)-1,cex=1.2)
  } else {
    if (md_log) { md_title <- paste(md_title,"(log scale)") } 
    temp_x <- seq(from=par("usr")[1] + (par("usr")[2] - par("usr")[1]) * .15,
                  to=par("usr")[2] - (par("usr")[2] - par("usr")[1]) * .15,
                  length.out=101)
    for (i in seq_along(idcol)) {
      rect(xleft=temp_x[i],xright=temp_x[i+1],
           ybottom=par("usr")[4] + (par("usr")[4] - par("usr")[3]) * .001,
           ytop=par("usr")[4] + strheight(md_title),
           col=idcol[i],border=NA,xpd=NA)
    }
    mtext(round(min(md),2),side=3,line=0,at=temp_x[1],adj=1.1)
    mtext(round(max(md),2),side=3,line=0,at=temp_x[101],adj=-0.1)
    mtext(md_title,side=3,line=1,at=temp_x[51],adj=.5,font=2,cex=1.2)
  }
}


# Metadata plots ----------
# ^ mdCompare -----------
plot_mdScatter <- function(MD,sel_clust,md_log) {
  temp_par <- par(no.readonly=T)
  layout(matrix(c(2,1,0,3),2),c(5,1),c(1,5))
  par(mar=c(3,3,0,0),mgp=2:0,cex=1.1)
  plot(MD[!MD$sel_cells,1:2],log=md_log,xlim=range(MD[,1]),ylim=range(MD[,2]),
       pch=21,col=alpha("black",0.2),bg=alpha("black",0.1))
  points(MD[MD$sel_cells,1:2],
         pch=21,col=alpha("red",0.4),bg=alpha("red",0.2))
  par(mar=c(0,3,1,0))
  boxplot(tapply(MD[,1],MD$sel_cells,c),log=sub("y","",md_log),
          horizontal=T,xaxt="n",yaxt="n",border=c("black","red"))
  if (any(MD$sel_cells)) {
    legend(x=switch(sub("y","",md_log),"x"=10^par("usr")[1],par("usr")[1]),y=par("usr")[4],
           xjust=0,yjust=0.2,xpd=NA,bty="n",pch=21,col="red",pt.bg=alpha("red",0.5),
           legend=paste("Cluster",names(sel_clust),sel_clust))
  }
  par(mar=c(3,0,0,1))
  boxplot(tapply(MD[,2],MD$sel_cells,c),log=sub("x","",md_log),
          horizontal=F,xaxt="n",yaxt="n",border=c("black","red"))
  par(temp_par)
}

plot_mdBoxplotX <- function(MD,sel_clust,md_log) {
  temp_par <- par(no.readonly=T)
  par(mar=c(3,3,2,1),mgp=2:0,cex=1.1)
  if (any(MD$sel_cells)) {
    temp1 <- tapply(MD[!MD$sel_cells,2],as.factor(MD[!MD$sel_cells,1]),c)
    temp2 <- tapply(MD[MD$sel_cells,2],as.factor(MD[MD$sel_cells,1]),c)
    plot(x=NULL,y=NULL,ylim=range(MD[,2]),
         xlim=c(0,length(levels(as.factor(MD[,1]))) * 3),
         log=md_log,xaxt="n",
         xlab=names(MD)[1],ylab=names(MD)[2])
    boxplot(temp1,add=T,xaxt="n",
            at=seq(1,length(levels(as.factor(MD[,1]))) * 3,by=3))
    boxplot(temp2,add=T,xaxt="n",border="red",
            at=seq(2,length(levels(as.factor(MD[,1]))) * 3,by=3))
    axis(side=1,at=seq(1.5,length(levels(as.factor(MD[,1]))) * 3,by=3),
         labels=names(temp1))
    legend(x=par("usr")[1],y=switch(sub("x","",md_log),"y"=10^par("usr")[4],par("usr")[4]),
           xjust=0,yjust=0.2,xpd=NA,bty="n",pch=0,col="red",pt.bg=alpha("red",0.5),
           legend=paste("Cluster",names(sel_clust),sel_clust))
    
  } else {
    boxplot(tapply(MD[,2],as.factor(MD[,1]),c),log=md_log,
            xlab=names(MD)[1],ylab=names(MD)[2])
  }
  par(temp_par)
}

plot_mdBoxplotY <- function(MD,sel_clust,md_log) {
  temp_par <- par(no.readonly=T)
  par(mar=c(3,3,2,1),mgp=2:0,cex=1.1)
  if (any(MD$sel_cells)) {
    temp1 <- tapply(MD[!MD$sel_cells,1],as.factor(MD[!MD$sel_cells,2]),c)
    temp2 <- tapply(MD[MD$sel_cells,1],as.factor(MD[MD$sel_cells,2]),c)
    plot(x=NULL,y=NULL,xlim=range(MD[,1]),
         ylim=c(0,length(levels(as.factor(MD[,2]))) * 3),
         log=md_log,yaxt="n",
         xlab=names(MD)[1],ylab=names(MD)[2])
    boxplot(temp1,add=T,yaxt="n",horizontal=T,
            at=seq(1,length(levels(as.factor(MD[,2]))) * 3,by=3))
    boxplot(temp2,add=T,yaxt="n",border="red",horizontal=T,
            at=seq(2,length(levels(as.factor(MD[,2]))) * 3,by=3))
    axis(side=2,at=seq(1.5,length(levels(as.factor(MD[,2]))) * 3,by=3),
         labels=names(temp1))
    legend(x=switch(sub("y","",md_log),"x"=10^par("usr")[1],par("usr")[1]),y=par("usr")[4],
           xjust=0,yjust=0.2,xpd=NA,bty="n",pch=0,col="red",pt.bg=alpha("red",0.5),
           legend=paste("Cluster",names(sel_clust),sel_clust))
    
  } else {
    boxplot(tapply(MD[,1],as.factor(MD[,2]),c),log=md_log,
            horizontal=T,xlab=names(MD)[1],ylab=names(MD)[2])
  }
  par(temp_par)
}


plot_mdCompare <- function(MD,mdX,mdY,sel_cells,sel_clust,md_log) {
  MD <- data.frame(MD[,c(mdX,mdY)])
  MD$sel_cells <- rownames(MD) %in% sel_cells
  if ("x" %in% md_log & !(is.factor(MD[,1]) | is.character(MD[,1]))) {
    tempLX <- "x"
    if (any(MD[,1] <= 0)) {
      names(MD)[1] <- paste(names(MD)[1],
                            paste0("(log scale: ",sum(MD[,1] <= 0),
                                   " values <= 0 omitted)"))
      MD <- MD[MD[,1] > 0,]
    } else {
      names(MD)[1] <- paste(names(MD)[1],"(log scale)")
    }
  } else {
    tempLX <- ""
  }
  if ("y" %in% md_log & !(is.factor(MD[,2]) | is.character(MD[,2]))) { 
    tempLY <- "y" 
    if (any(MD[,2] <= 0)) {
      names(MD)[2] <- paste(names(MD)[2],
                            paste0("(log scale: ",sum(MD[,2] <= 0),
                                   " values <= 0 omitted)"))
      MD <- MD[MD[,2] > 0,]
    } else {
      names(MD)[2] <- paste(names(MD)[2],"(log scale)")
    } 
  } else {
    tempLY <- ""
  }
  md_log <- paste(c(tempLX,tempLY),collapse="")
  
  if ((is.factor(MD[,1]) | is.character(MD[,1])) &
      (is.factor(MD[,2]) | is.character(MD[,2]))) {
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,"This figure is not designed to compare to categorical variables.")
  } else if (is.factor(MD[,1]) | is.character(MD[,1])) {
    plot_mdBoxplotX(MD,sel_clust,md_log)
  } else if (is.factor(MD[,2]) | is.character(MD[,2])) {
    plot_mdBoxplotY(MD,sel_clust,md_log)
  } else {
    plot_mdScatter(MD,sel_clust,md_log)
  }
  
}


# plot_mdCompare(MD=getMD(inD),
#                mdX="total_counts",
#                mdY="cyclonePhases",
#                sel_cells=Clusters(sCVdL$res.1.4) == "1", # rep(F,nrow(getMD(inD))),
#                sel_clust=attr(Clusters(sCVdL$res.1.4),"ClusterNames")[1],
#                md_log=c("x","y"))


# ^ mdPerClust -------
plot_mdBarplot <- function(MD,opt) {
  temp_par <- par(no.readonly=T)
  id0 <- as.factor(MD[,1])
  id <- switch(opt,
               "relative"=tapply(id0,MD$cl,function(X) table(X) / length(X)),
               "absolute"=tapply(id0,MD$cl,table))
  if (is.list(id)) { id <- do.call(cbind,id) }
  idylab <- switch(opt,
                   "relative"="Proportion of cells per cluster",
                   "absolute"="Number of cells per cluster")
  if (length(levels(id0)) <= 8) {
    idcol <- RColorBrewer::brewer.pal(length(levels(id0)),"Dark2")[1:length(levels(id0))]
  } else {
    idcol <- rainbow2(length(levels(id0)))
  }
  par(mar=c(3,3,ceiling(length(levels(id0))/4)+.5,1),mgp=2:0,cex=1.1)
  barplot(id,col=idcol,ylab=idylab,xlab="Clusters",yaxt="n",mgp=c(2,0,0),
          legend.text=levels(id0),font=2,
          args.legend=list(x=par("usr")[2],y=par("usr")[4],
                           xjust=1,yjust=0.2,xpd=NA,ncol=4,bty="n"))
  axis(2)
  barplot(rep(par("usr")[3]*4.6,length(levels(MD$cl))),add=T,
          col=alpha(rainbow2(length(levels(MD$cl))),0.5),
          border=alpha(rainbow2(length(levels(MD$cl))),0.5))
  abline(h=0)
  mtext(names(MD)[1],side=3,adj=0,font=2,line=ceiling(length(levels(id0))/4)-.8,cex=1.2)
  par(temp_par)
}

plot_mdBoxplot <- function(MD,opt) {
  temp_par <- par(no.readonly=T)
  par(mar=c(3,3,2,1),mgp=2:0,cex=1.1)
  boxplot(tapply(MD[,1],MD$cl,c),log=opt,
          ylab=names(MD)[1],xlab="Clusters",
          border=rainbow2(length(levels(MD$cl))),
          col=alpha(rainbow2(length(levels(MD$cl))),0.3))
  par(temp_par)
}

plot_mdPerClust <- function(MD,sel,cl,opt) {
  MD <- MD[sel]
  MD$cl <- cl
  if ("y" %in% opt & !(is.factor(MD[,1]) | is.character(MD[,1]))) { 
    if (any(MD[,1] <= 0)) {
      names(MD)[1] <- paste(names(MD)[1],
                            paste0("(log scale: ",sum(MD[,1] <= 0),
                                   " values <= 0 omitted)"))
      MD <- MD[MD[,1] > 0,]
    } else {
      names(MD)[1] <- paste(names(MD)[1],"(log scale)")
    } 
  }
  if (is.factor(MD[,1]) | is.character(MD[,1])) {
    plot_mdBarplot(MD,opt)
  } else {
    plot_mdBoxplot(MD,opt)
  }
}

# plot_mdPerClust(MD=getMD(inD),sel="cyclonePhases",cl=Clusters(sCVdL$res.0.8),opt="absolute")


# DE gene dotplot -----------
dotplotDEgenes <- function(sCVd,DEtype,FDRthresh) {
  if (missing(FDRthresh)) { FDRthresh <- 1 }
  if (DEtype == "DEvsRest") {
    return(lapply(DEvsRest(sCVd),function(X) {
      temp <- X[which(X$FDR <= FDRthresh),"FDR",drop=F]
      out <- unlist(temp,use.names=F)
      names(out) <- rownames(temp)
      return(sort(out))
    }))
  } else if (DEtype == "DEneighb") {
    outL <- lapply(DEneighb(sCVd,FDRthresh), function(X) {
      if (nrow(X) < 1) { return(numeric(0)) }
      out <- X[,grep("^FDR_",names(X))]
      names(out) <- rownames(X)
      return(sort(out))
    })
    names(outL) <- levels(Clusters(sCVd))
    return(outL)
  } else if (DEtype == "DEmarker") {
    outL <- lapply(DEmarker(sCVd,FDRthresh), function(X) {
      if (nrow(X) < 1) { return(numeric(0)) }
      out <- apply(X[,grep("^FDR_",names(X)),drop=F],1,max)
      return(sort(out))
    })
    return(outL)
  }
}
# returns A named list of named numeric vectors, one entry for each cluster.
# Names are genes passing FDR threshold, and value is (maximum) FDR of comparison.


plot_deDotplot <- function(sCVd,DEgenes,DEnum) {
  # ^ Setup ----
  heatGenes <- unique(unlist(lapply(DEgenes,function(X) names(X)[1:DEnum])))
  heatGenes <- heatGenes[!is.na(heatGenes)]
  
  if (is.null(heatGenes)) {
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,
         "No genes were statistically significant at the current false discovery rate")
    return(invisible())
  }
  
  temp_DR <- sapply(ClustGeneStats(sCVd),function(X) X[heatGenes,"DR"])
  if (is.vector(temp_DR)) {
    temp_DR <- matrix(temp_DR,1,dimnames=list(NULL,names(temp_DR))) 
  }
  temp_MDGE <- sapply(ClustGeneStats(sCVd),function(X) X[heatGenes,"MDGE"])
  if (is.vector(temp_MDGE)) {
    temp_MDGE <- matrix(temp_MDGE,1,dimnames=list(NULL,names(temp_MDGE))) 
  }
  rownames(temp_DR) <- rownames(temp_MDGE) <- heatGenes
  
  if (nrow(temp_DR) > 1) {
    hG <- hclust(dist(temp_DR),"complete")
  } else {
    hG <- list(order=1)
  }
  
  #  if (length(ClustGeneStats(sCVd)) > 2) {
  hC <- hclust(as.dist(DEdist(sCVd)),"single")
  # } else {
  #    hC <- hclust(dist(t(temp_DR)))
  #  }  
  clustCols <- rainbow2(length(levels(Clusters(sCVd))))
  
  dC <- dendrapply(as.dendrogram(hC),function(X) {
    if (is.leaf(X)) {
      attr(X,"edgePar") <- list(
        lwd=2,
        col=clustCols[which(attr(X,"label") == levels(Clusters(sCVd)))]
      )
      attr(X,"nodePar") <- list(
        pch=NA,lab.font=2,lab.cex=1.2,
        lab.col=clustCols[which(attr(X,"label") == levels(Clusters(sCVd)))])
      if (attr(X,"label") != "Unselected") {
        if (attr(X,"label") %in% names(DEgenes)) {
          attr(X,"label") <- paste0(attr(X,"label"),": ",
                                    length(DEgenes[[attr(X,"label")]])," DE")
        } else {
          attr(X,"label") <- paste0(
            attr(X,"label"),": ",
            length(DEgenes[[which(attr(X,"label") ==
                                    sapply(strsplit(names(DEgenes),"-"),
                                           function(X) X[1]))]]),
            " DE")
        }
      }
    }
    return(X)
  })
  
  tempLabCol <- ClustGeneStats(sCVd)[[1]][heatGenes,"genes"]
  DR <- temp_DR[hG$order,hC$order,drop=F]
  temp <- range(sapply(ClustGeneStats(sCVd),function(X) X[,"MDGE"]))
  temp <- seq(temp[1],temp[2],length.out=101)
  MDGE <- findInterval(as.vector(temp_MDGE[hG$order,hC$order]),
                       vec=temp,all.inside=T)
  
  # ^ Plot dotplot ----
  temp_par <- par(no.readonly=T)
  if (length(levels(Clusters(sCVd))) <= 1) {
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,paste("Heatmap cannot be computed",
                     "with less than two clusters.",sep="\n"))
  } else if (length(heatGenes) < 1) {
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,paste("There are no differentially expressed genes at",
                     "false discovery rate threshold."))
  } else {
    layout(matrix(c(0,2,3,1),2),widths=c(1,5),heights=c(1,5))
    par(mar=c(9,0,0,.5))
    plot(x=NULL,y=NULL,xlim=c(0.5,nrow(DR)+.5),ylim=c(0.5,ncol(DR)+.5),
         xaxs="i",yaxs="i",xaxt="n",yaxt="n",xlab=NA,ylab=NA,bty="n")
    abline(v=1:nrow(DR),col="grey90")
    symbols(x=rep(1:nrow(DR),ncol(DR)),
            y=as.vector(sapply(1:ncol(DR),function(X) rep(X,nrow(DR)))),
            circles=as.vector(DR)/2,inches=F,add=T,xpd=NA,
            fg=viridis::viridis(100,d=-1)[MDGE],
            bg=viridis::viridis(100,d=-1)[MDGE])
    axis(side=1,at=1:nrow(DR),lwd=0,labels=tempLabCol[hG$order],las=2,cex.axis=1.2)
    
    # Legend:
    tx0 <- par("usr")[1]
    tx <- (par("usr")[2] - par("usr")[1])
    ty0 <- par("usr")[3]
    ty <- par("usr")[4] - par("usr")[3]
    segments(x0=tx0 - seq(.15,.03,length.out=1000) * tx,
             y0=ty0 - 0.02 * ty,y1=ty0 - 0.05 * ty,
             col=viridis::viridis(1000,d=-1),xpd=NA)
    text(x=tx0 - c(.15,.09,.03) * tx,
         y=ty0 - c(0.035,0.02,0.035) * ty,
         labels=c(round(min(temp_MDGE),2),
                  "Mean detected expression",
                  round(max(temp_MDGE),2)),pos=2:4,xpd=NA)
    symbols(x=tx0 - c(.15,.09,.03) * tx,
            y=ty0 - rep(.14,3) * ty,add=T,xpd=NA,
            circles=c(0.25,0.5,0.75)/2,inches=F,bg="black")
    text(x=tx0 - c(.149,.089,0.029,.09) * tx,
         y=ty0 - c(rep(.23,3),.26) * ty,xpd=NA,
         labels=c("25%","50%","75%","Detection Rate"))
    
    
    par(mar=c(9,0,0,7))
    plot(dC,horiz=T,xpd=NA,
         ylim=c(0.5,length(hC$order)+.5),yaxs="i",yaxt="n")
    
    par(mar=c(0,0,0,.5))
    if (class(hG) == "hclust") {
      plot(as.dendrogram(hG),leaflab="none",
           xlim=c(0.5,length(hG$order)+.5),xaxs="i",yaxt="n")
    }
  }
  par(temp_par)
}

#plot_deDotplot(sCVdL$res.0.4,NULL,"deMarker",5,0.001)


# Cluster-wise gene expression scatterplot -----
# ^ dot types
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
doubleDot <- function(col1,col2) {
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



plot_clusterGenes_markers <- function(sCVd,selClust,cellMarkersS,cellMarkersU) {
  cellMarkCols <- rainbow2(length(cellMarkersU))
  par(mar=c(3,3,3,20),mgp=2:0)
  if (selClust == "") {
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,paste("Click a cell from a cluster on the tSNE plot above",
                     "or select a cluster from the drop-down list above left",
                     "to see gene expression for that cluster.",sep="\n"))
  } else {
    CGS <- ClustGeneStats(sCVd)[[selClust]]
    temp_ylab <- switch(as.character(Param(sCVd,"exponent") == exp(1)),
                        "TRUE"="(natural log scale)",
                        "FALSE"=paste0("(log",Param(sCVd,"exponent")," scale)"))
    plot(MDGE~DR,
         data=CGS[!((CGS$cMu | CGS$cMs) & CGS$overCut),],
         xlim=range(CGS$DR),ylim=range(CGS$MDGE),
         col=alpha("black",0.2),pch=20,
         xlab="Proportion of cells in which gene was detected",
         ylab=paste("Mean normalized gene expression where detected",temp_ylab))
    title(paste0("Cluster ", selClust,": ",attr(Clusters(sCVd),"ClusterNames")[selClust]),cex=1.2)
    mtext(paste("Cells:",sum(Clusters(sCVd) == selClust),
                "   Genes detected:",sum(CGS$DR > 0)),side=3,line=0,cex=0.9)
    box(col=rainbow2(length(levels(Clusters(sCVd))))[selClust],lwd=2)
    for (x in which(CGS$cMu)) {
      TeachingDemos::my.symbols(x=CGS$DR[x],y=CGS$MDGE[x],
                                symb=singleDot,inches=0.1,
                                MoreArgs=list(
                                  col1=cellMarkCols[which(
                                    sapply(cellMarkersU,function(X) 
                                      CGS$genes[x] %in% X)
                                  )]
                                ))
    }
    for (x in which(CGS$cMs)) {
      temp <- unlist(strsplit(names(which(sapply(cellMarkersS,function(X) 
        CGS$genes[x] %in% X))),"&"))
      TeachingDemos::my.symbols(x=CGS$DR[x],
                                y=CGS$MDGE[x],
                                symb=doubleDot,inches=0.1,
                                MoreArgs=list(col1=cellMarkCols[as.integer(temp[1])],
                                              col2=cellMarkCols[as.integer(temp[2])]))
    }
    tempLabels <- spreadLabels(CGS[(CGS$cMu | CGS$cMs) & CGS$overCut,"DR"],
                               CGS[(CGS$cMu | CGS$cMs) & CGS$overCut,"MDGE"],
                               CGS[(CGS$cMu | CGS$cMs) & CGS$overCut,"genes"],
                               str.cex=1.2,str.font=2)
    rownames(tempLabels) <- CGS[(CGS$cMu | CGS$cMs) & CGS$overCut,"genes"]                        
    for (gn in CGS[CGS$cMu & CGS$overCut,"genes"]) {
      rect(xleft=tempLabels[gn,1] - strwidth(gn,cex=1.2,font=2) * .5,
           xright=tempLabels[gn,1] + strwidth(gn,cex=1.2,font=2) * .5,
           ybottom=tempLabels[gn,2] - strheight(gn,cex=1.2,font=2) * .5,
           ytop=tempLabels[gn,2] + strheight(gn,cex=1.2,font=2) * .5,
           border=NA,col=alpha("white",0.5))
      text(tempLabels[gn,,drop=F],labels=gn,cex=1.2,font=2,
           col=cellMarkCols[which(sapply(cellMarkersU,function(X) 
             gn %in% X))])
    }
    for (gn in CGS[CGS$cMs & CGS$overCut,"genes"]) {
      rect(xleft=tempLabels[gn,1] - strwidth(gn,cex=1.2,font=2) * .5,
           xright=tempLabels[gn,1] + strwidth(gn,cex=1.2,font=2) * .5,
           ybottom=tempLabels[gn,2] - strheight(gn,cex=1.2,font=2) * .5,
           ytop=tempLabels[gn,2] + strheight(gn,cex=1.2,font=2) * .5,
           border=NA,col=alpha("white",0.5))
      text(tempLabels[gn,,drop=F],labels=gn,cex=1.2,font=2,
           col=cellMarkCols[as.integer(temp[2])])
    }
    legend(x=1.05,y=max(CGS$MDGE),xpd=NA,bty="n",ncol=1,
           pch=19,col=cellMarkCols,legend=names(cellMarkersU))
  }
}


plot_clusterGenes_DEgenes <- function(sCVd,selClust,DEgenes,DEnum,DEtype) {
  par(mar=c(3,3,3,20),mgp=2:0)
  if (selClust == "") {
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,paste("Click a cell from a cluster on the tSNE plot above",
                     "or select a cluster from the drop-down list above left",
                     "to see gene expression for that cluster.",sep="\n"))
  } else {
    CGS <- ClustGeneStats(sCVd)[[selClust]]
    temp_ylab <- switch(as.character(Param(sCVd,"exponent") == exp(1)),
                        "TRUE"="(natural log scale)",
                        "FALSE"=paste0("(log",Param(sCVd,"exponent")," scale)"))
    plot(MDGE~DR,
         data=CGS[!rownames(CGS) %in% names(DEgenes[[selClust]])[1:DEnum],],
         xlim=range(CGS$DR),ylim=range(CGS$MDGE),
         col=alpha("black",0.2),pch=20,
         xlab="Proportion of cells in which gene was detected",
         ylab=paste("Mean normalized gene expression where detected",temp_ylab))
    title(paste0("Cluster ", selClust,": ",attr(Clusters(sCVd),"ClusterNames")[selClust]),cex=1.2)
    mtext(paste("Cells:",sum(Clusters(sCVd) == selClust),
                "   Genes detected:",sum(CGS$DR > 0)),side=3,line=0,cex=0.9)
    box(col=rainbow2(length(levels(Clusters(sCVd))))[selClust],lwd=2)
    if (length(DEgenes[[selClust]]) > 0) {
      DEG <- names(DEgenes[[selClust]])[1:DEnum]
      DEG <- DEG[!is.na(DEG)]
      points(x=CGS[DEG,"DR"],y=CGS[DEG,"MDGE"],
             pch=16,cex=1.2,col="firebrick2")
      if (any(CGS[DEG,"overCut"])) {
        tempLabels <- spreadLabels(x=CGS[DEG,"DR"][CGS[DEG,"overCut"]],
                                   y=CGS[DEG,"MDGE"][CGS[DEG,"overCut"]],
                                   label=CGS[DEG,"genes"][CGS[DEG,"overCut"]],
                                   str.cex=1.2,str.font=2)
        rect(xleft=tempLabels[,1] - 
               strwidth(CGS[DEG,"genes"][CGS[DEG,"overCut"]],cex=1.2,font=2) * .5,
             xright=tempLabels[,1] + 
               strwidth(CGS[DEG,"genes"][CGS[DEG,"overCut"]],cex=1.2,font=2) * .5,
             ybottom=tempLabels[,2] - 
               strheight(CGS[DEG,"genes"][CGS[DEG,"overCut"]],cex=1.2,font=2) * .5,
             ytop=tempLabels[,2] +
               strheight(CGS[DEG,"genes"][CGS[DEG,"overCut"]],cex=1.2,font=2) * .5,
             border=NA,col=alpha("white",0.5))
        text(tempLabels,cex=1.2,font=2,col="firebrick2",
             labels=CGS[DEG,"genes"][CGS[DEG,"overCut"]])
      }
    }
    temp_n <- length(DEgenes[[selClust]])
    temp_lab <- switch(DEtype,
                       DEvsRest=" DE genes vs rest of cells in sample",
                       DEmarker=" marker genes",
                       deNeighb=" DE genes vs nearest neighbouring cluster")
    legend("top",bty="n",pch=16,col="firebrick2",
           legend=paste0(temp_n,temp_lab," (showing top ",
                         min(temp_n,DEnum),")"))
    
  }
}

plot_clusterGenes_search <- function(sCVd,selClust,GOI) {
  par(mar=c(3,3,3,20),mgp=2:0)
  if (selClust == "") {
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,paste("Click a cell from a cluster on the tSNE plot above",
                     "or select a cluster from the drop-down list above left",
                     "to see gene expression for that cluster.",sep="\n"))
  } else {
    CGS <- ClustGeneStats(sCVd)[[selClust]]
    temp_ylab <- switch(as.character(Param(sCVd,"exponent") == exp(1)),
                        "TRUE"="(natural log scale)",
                        "FALSE"=paste0("(log",Param(sCVd,"exponent")," scale)"))
    plot(MDGE~DR,
         data=CGS[!rownames(CGS) %in% GOI,],
         col=alpha("black",0.2),pch=20,
         xlim=range(CGS$DR),ylim=range(CGS$MDGE),
         xlab="Proportion of cells in which gene was detected",
         ylab=paste("Mean normalized gene expression where detected",temp_ylab))
    title(paste0("Cluster ", selClust,": ",attr(Clusters(sCVd),"ClusterNames")[selClust]),cex=1.2)
    mtext(paste("Cells:",sum(Clusters(sCVd) == selClust),
                "   Genes detected:",sum(CGS$DR > 0)),side=3,line=0,cex=0.9)
    box(col=rainbow2(length(levels(Clusters(sCVd))))[selClust],lwd=2)
    if (length(GOI) > 0) {
      points(x=CGS[GOI,"DR"],y=CGS[GOI,"MDGE"],
             pch=16,cex=1.2,col="firebrick2")
      if (any(CGS[GOI,"overCut"])) {
        tempLabels <- spreadLabels(x=CGS[GOI,"DR"][CGS[GOI,"overCut"]],
                                   y=CGS[GOI,"MDGE"][CGS[GOI,"overCut"]],
                                   label=CGS[GOI,"genes"][CGS[GOI,"overCut"]],
                                   str.cex=1.2,str.font=2)
        rect(xleft=tempLabels[,1] - 
               strwidth(CGS[GOI,"genes"][CGS[GOI,"overCut"]],cex=1.2,font=2) * .5,
             xright=tempLabels[,1] + 
               strwidth(CGS[GOI,"genes"][CGS[GOI,"overCut"]],cex=1.2,font=2) * .5,
             ybottom=tempLabels[,2] - 
               strheight(CGS[GOI,"genes"][CGS[GOI,"overCut"]],cex=1.2,font=2) * .5,
             ytop=tempLabels[,2] +
               strheight(CGS[GOI,"genes"][CGS[GOI,"overCut"]],cex=1.2,font=2) * .5,
             border=NA,col=alpha("white",0.5))
        text(tempLabels,cex=1.2,font=2,col="firebrick2",
             labels=CGS[GOI,"genes"][CGS[GOI,"overCut"]])
      }
    }
  }
}

# Gene search function -------
geneSearch <- function(txt,st,CGS) {
  if (length(txt) < 1) { txt <- ""}
  geneNames <- rownames(CGS)
  names(geneNames) <- toupper(CGS$genes)
  temp <- switch(st,
                 comma={
                   temp_in <- strsplit(txt,split="[\\s,]",perl=T)[[1]]
                   temp_out <- geneNames[toupper(temp_in)]
                   names(temp_out) <- CGS[temp_out,"genes"]
                   return(temp_out)
                 },
                 regex={
                   temp_in <- grep(txt,names(geneNames),ignore.case=T)
                   temp_out <- geneNames[temp_in]
                   names(temp_out) <- CGS[temp_out,"genes"]
                   return(temp_out)
                 })
  if (length(temp) > 0) {
    return(temp)
  } else {
    return(switch(st,
                  comma={
                    temp_in <- strsplit(txt,split="[\\s,]",perl=T)[[1]]
                    return(geneNames[which(toupper(geneNames) %in% toupper(temp_in))])
                  },
                  regex=grep(txt,geneNames,value=T,ignore.case=T)))
  }
}

# Gene expression boxplots --------
plot_GEboxplot <- function(nge,sCVd,gene,geneName,opts) {
  if (gene == "") {
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,paste("Select a gene by either clicking on the plot above",
                     "or searching for genes of interest in the search bar above,",
                     "then pick the gene from the list just above this figure",
                     "to see a comparison of that gene's expression across all clusters.",
                     sep="\n"))
  } else {
    # ^ setup -----
    hC <- hclust(as.dist(DEdist(sCVd)),"single")
    temp_ylab <- switch(as.character(Param(sCVd,"exponent") == exp(1)),
                        "TRUE"="(natural log scale)",
                        "FALSE"=paste0("(log",Param(sCVd,"exponent")," scale)"))
    temp_pos <- switch(as.character(length(levels(Clusters(sCVd))) > 1),
                       "TRUE"=hC$order,"FALSE"=1)
    if ("sct" %in% opts) {
      bxpCol <- rainbow2(length(levels(Clusters(sCVd))),.2)
    } else {
      bxpCol <- rainbow2(length(levels(Clusters(sCVd))),.8)
    }
    
    # ^ plot boxplot -----
    temp_par <- par(no.readonly=T)
    layout(matrix(2:1,nrow=2),heights=c(1,4))
    par(mar=c(3,3,0,3),mgp=2:0)
    suppressWarnings(boxplot(
      vector("list",length(levels(Clusters(sCVd))[levels(Clusters(sCVd)) != "Unselected"])),
      ylim=range(nge[gene,]),xaxt="n",xlab=NA,
      ylab=paste(gene,"normalized gene expression",temp_ylab)
    ))
    mtext(levels(Clusters(sCVd))[temp_pos],side=1,line=0,at=seq_along(temp_pos))
    mtext("Clusters, ordered by heatmap dendrogram",side=1,line=1)
    if (is.null(geneName)) { 
      mtext(paste(gene,collapse="\n"),side=1,line=2,font=2) 
    } else {
      mtext(paste(paste0(names(geneName),": ",geneName),collapse="\n"),
            side=1,line=2,font=2) 
    }
    for (i in temp_pos) {
      boxplot(nge[gene,Clusters(sCVd) %in% levels(Clusters(sCVd))[i]],add=T,
              at=which(temp_pos == i),col=bxpCol[i],outline=F)
      if ("sct" %in% opts) {
        points(jitter(rep(which(temp_pos == i),
                          sum(Clusters(sCVd) %in% levels(Clusters(sCVd))[i])),
                      amount=.2),
               nge[gene,Clusters(sCVd) %in% levels(Clusters(sCVd))[i]],
               pch=20,col=rainbow2(length(levels(Clusters(sCVd))),.4)[i])
      }
    }
    if ("dr" %in% opts) {
      points(x=seq_along(ClustGeneStats(sCVd)),
             y=sapply(ClustGeneStats(sCVd)[temp_pos],function(X) X[gene,"DR"]) * 
               max(nge[gene,]) + min(nge[gene,]),
             pch="-",cex=2)
      axis(side=4,at=seq(0,1,.25) * max(nge[gene,]) + min(nge[gene,]),
           labels=paste0(seq(0,1,.25) * 100,"%"))
      mtext(side=4,line=2,text="- Gene detection rate per cluster")
    }
    if (length(temp_pos) > 1) { 
      par(new=F,mar=c(0,3,1,3))
      plot(as.dendrogram(hC),leaflab="none") 
    }
    par(temp_par)
  }
}


# Cluster comparisons --------
compareClusts_DF <- function(sCVd,clA,clB,dataType) {
  if (dataType %in% c("MGE","MDGE","DR")) {
    loc1 <- c(paste(clA,clB,sep="-"),paste(clB,clA,sep="-"))
    loc <- loc1[loc1 %in% names(DEcombn(sCVd))]
    loc1 <- which(loc1 %in% names(DEcombn(sCVd)))
    if (loc1 == 2) { loc1 <- -1 }
    tempW <- DEcombn(sCVd)[[loc]]$Wstat - 
      DEcombn(sCVd)[[loc]]$Wstat[which.max(DEcombn(sCVd)[[loc]]$pVal)]
    temp <- data.frame(x_diff=ClustGeneStats(sCVd)[[clA]][,dataType] - 
                         ClustGeneStats(sCVd)[[clB]][,dataType],
                       y_mean=rowMeans(cbind(ClustGeneStats(sCVd)[[clA]][,dataType],
                                             ClustGeneStats(sCVd)[[clB]][,dataType])),
                       logGER=DEcombn(sCVd)[[loc]]$logGER,
                       FDR=DEcombn(sCVd)[[loc]]$FDR,
                       dir=c(clB,clA)[(tempW * loc1 > 0) + 1])
    rownames(temp) <- rownames(ClustGeneStats(sCVd)[[clA]])
    return(temp)
  } else if (dataType %in% c("GERvDDR","logGER","dDR")) {
    loc1 <- which(c(paste(clA,clB,sep="-"),paste(clB,clA,sep="-")) %in% names(DEcombn(sCVd)))
    if (loc1 == 2) { loc1 <- -1 }
    loc <- which(names(DEcombn(sCVd)) %in% c(paste(clA,clB,sep="-"),paste(clB,clA,sep="-")))
    temp <- DEcombn(sCVd)[[loc]][,c("logGER","dDR","FDR")]
    temp <- as.data.frame(mapply("*",temp,c(loc1,loc1,1))) 
    rownames(temp) <- rownames(ClustGeneStats(sCVd)[[clA]])
    tempW <- DEcombn(sCVd)[[loc]]$Wstat - 
      DEcombn(sCVd)[[loc]]$Wstat[which.max(DEcombn(sCVd)[[loc]]$pVal)]
    temp$dir <- c(clB,clA)[(tempW * loc1 > 0) + 1]
    return(temp)
  } 
}

plot_compareClusts_MAplot <- function(sCVd,clA,clB,dataType,labType,labNum,labGenes,labAngle) {
  # ^ setup -----
  CGS <- compareClusts_DF(sCVd,clA,clB,dataType)
  temp_exp <- switch(as.character(Param(sCVd,"exponent") == exp(1)),
                     "TRUE"="(natural log scale)",
                     "FALSE"=paste0("(log",Param(sCVd,"exponent")," scale)"))
  temp_label <- switch(dataType,
                       "MGE"=paste("mean normalized gene expression",temp_exp),
                       "MDGE"=paste("mean normalized gene expression where detected",temp_exp),
                       "DR"="proportion of cells in which gene was detected")
  if ("flip" %in% labAngle) {
    temp_adjB <- c(-0.1,0.5)
    temp_adjA <- c(1.1,0.5)
  } else {
    temp_adjB <- c(1.1,0.5)
    temp_adjA <- c(-0.1,0.5)
  }
  if (labType == "diff") {
    gnA <- rownames(head(CGS[order(CGS$x_diff,decreasing=T),],labNum))
    gnB <- rownames(tail(CGS[order(CGS$x_diff,decreasing=T),],labNum))
  } else if (labType == "de") {
    ts <- order(CGS$FDR,na.last=T)
    gnA <- rownames(CGS)[ts[CGS[ts,"dir"] == clA][1:labNum]]
    gnB <- rownames(CGS)[ts[CGS[ts,"dir"] == clB][1:labNum]]
  }
  
  # ^ plot -----
  par(mar=c(3,3,3.5,1),mgp=c(2,1,0))
  plot(y_mean~x_diff,data=CGS,
       xlab=paste0("Difference in ",temp_label," (",clA," - ",clB,")"),
       ylab=paste0("Average of ",temp_label," between ",clA," & ",clB),
       main=paste0("Modified MA plot of ",
                   switch(dataType,
                          "MGE"="mean gene expression",
                          "MDGE"="mean detected gene expression",
                          "DR"="detection rate"),
                   " (",clA," vs. ",clB,")"),
       pch=20,col=alpha("black",0.3))
  abline(v=0,col="gray50")
  lines(x=c(par("usr")[2],par("usr")[2]),y=c(par("usr")[3],par("usr")[4]),lwd=2,xpd=NA,
        col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clA)])
  lines(x=c(par("usr")[1],par("usr")[1]),y=c(par("usr")[3],par("usr")[4]),xpd=NA,lwd=2,
        col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clB)])
  if (labType == "search") {
    if (length(labGenes) > 0) {
      points(y_mean~x_diff,data=CGS[labGenes,],pch=16,col=alpha("firebrick2",0.8))
      tempLabel <- spreadLabels(x=CGS[labGenes,"x_diff"],y=CGS[labGenes,"y_mean"],
                                label=labGenes,str.cex=1.2,str.font=2)
      text(tempLabel,labels=labGenes,col="firebrick2",cex=1.2,font=2)
    }
  } else {
    points(y_mean~x_diff,data=CGS[gnA,],pch=16,
           col=rainbow2(length(levels(Clusters(sCVd))),.8)[which(levels(Clusters(sCVd)) == clA)])
    points(y_mean~x_diff,data=CGS[gnB,],pch=16,
           col=rainbow2(length(levels(Clusters(sCVd))),.8)[which(levels(Clusters(sCVd)) == clB)])
    tempLabel <- spreadLabels(x=CGS[c(gnA,gnB),"x_diff"],y=CGS[c(gnA,gnB),"y_mean"],
                              label=c(gnA,gnB),str.cex=1.2,str.font=2)
    rownames(tempLabel) <- c(gnA,gnB)
    text(tempLabel[gnA,],labels=gnA,cex=1.2,font=2,
         col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clA)])
    text(tempLabel[gnB,],labels=gnB,cex=1.2,font=2,
         col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clB)])
  }
  mtext(paste("Higher in",clA),side=1,line=-1.1,adj=.99,font=2,
        col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clA)])
  mtext(paste("Higher in",clB),side=1,line=-1.1,adj=0.01,font=2,
        col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clB)])
  mtext(paste(
    paste("Cosine similarity of comparison:",
          round(cosineSim(ClustGeneStats(sCVd)[[clA]][,dataType], 
                          ClustGeneStats(sCVd)[[clB]][,dataType]),2)),
    paste("Spearman's Rho of comparison:",
          round(cor(x=ClustGeneStats(sCVd)[[clA]][,dataType], 
                    y=ClustGeneStats(sCVd)[[clB]][,dataType],
                    method="spearman"),2)),
    sep="  -  "),
    side=3)
} 

plot_compareClusts_DEscatter <- function(sCVd,clA,clB,dataType,labType,
                                         labTypeDiff,labNum,labGenes,labAngle) {
  # ^ setup -----
  CGS <- compareClusts_DF(sCVd,clA,clB,dataType)
  temp_exp <- switch(as.character(Param(sCVd,"exponent") == exp(1)),
                     "TRUE"="(natural log scale)",
                     "FALSE"=paste0("(log",Param(sCVd,"exponent")," scale)"))
  if ("flip" %in% labAngle) {
    temp_adjB <- c(-0.1,0.5)
    temp_adjA <- c(1.1,0.5)
  } else {
    temp_adjB <- c(1.1,0.5)
    temp_adjA <- c(-0.1,0.5)
  }
  if (labType == "diff") {
    gnA <- rownames(head(CGS[order(CGS[[labTypeDiff]],decreasing=T),],labNum))
    gnB <- rownames(tail(CGS[order(CGS[[labTypeDiff]],decreasing=T),],labNum))
  } else if (labType == "de") {
    ts <- order(CGS$FDR,na.last=T)
    gnA <- rownames(CGS)[ts[CGS[ts,"dir"] == clA][1:labNum]]
    gnB <- rownames(CGS)[ts[CGS[ts,"dir"] == clB][1:labNum]]
  }
  # Adding a colour scale for FDR
  # CGS <- CGS[order(CGS$FDR,decreasing=T,na.last=F),]
  # temp_col <- viridis(100,alpha=0.3,direction=-1)[cut(-log10(CGS$FDR),100,labels=F)]
  # temp_col[is.na(temp_col)] <- alpha("grey90",0.3)
  
  # ^ plot -----
  par(mar=c(3,3,3.5,1),mgp=c(2,1,0))
  plot(logGER~dDR,data=CGS,
       xlab=paste0("Difference in detection rate (",clA," - ",clB,")"),
       ylab=paste0("Gene expression ratio (",clA," : ",clB,") ",temp_exp),
       main=paste0("Expression difference effect sizes (",clA," vs. ",clB,")"),
       pch=20,col=alpha("black",0.3)) # col=temp_col)
  abline(v=0,h=0,col="gray50")
  lines(x=c(par("usr")[2],par("usr")[2]),y=c(0,par("usr")[4]),lwd=2,xpd=NA,
        col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clA)])
  lines(x=c(0,par("usr")[2]),y=c(par("usr")[4],par("usr")[4]),lwd=2,xpd=NA,
        col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clA)])
  lines(x=c(par("usr")[1],par("usr")[1]),y=c(par("usr")[3],0),xpd=NA,lwd=2,
        col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clB)])
  lines(x=c(par("usr")[1],0),y=c(par("usr")[3],par("usr")[3]),xpd=NA,lwd=2,
        col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clB)])
  if (labType == "search") {
    if (length(labGenes) > 0) {
      points(logGER~dDR,data=CGS[labGenes,],pch=16,col=alpha("firebrick2",0.8))
      tempLabel <- spreadLabels(CGS[labGenes,"dDR"],CGS[labGenes,"logGER"],
                                label=labGenes,str.cex=1.2,str.font=2)
      text(tempLabel,labels=labGenes,col="firebrick2",cex=1.2,font=2)
    }
  } else {
    points(logGER~dDR,data=CGS[gnA,],pch=16,
           col=rainbow2(length(levels(Clusters(sCVd))),.8)[which(levels(Clusters(sCVd)) == clA)])
    points(logGER~dDR,data=CGS[gnB,],pch=16,
           col=rainbow2(length(levels(Clusters(sCVd))),.8)[which(levels(Clusters(sCVd)) == clB)])
    tempLabel <- spreadLabels(CGS[c(gnA,gnB),"dDR"],CGS[c(gnA,gnB),"logGER"],
                              label=c(gnA,gnB),str.cex=1.2,str.font=2)
    rownames(tempLabel) <- c(gnA,gnB)
    text(tempLabel[gnA,],labels=gnA,cex=1.2,font=2,
         col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clA)])
    text(tempLabel[gnB,],labels=gnB,cex=1.2,font=2,
         col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clB)])
  }
  mtext(paste("Higher in",clA),side=3,line=-1.1,adj=.99,font=2,
        col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clA)])
  mtext(paste("Higher in",clB),side=1,line=-1.1,adj=0.01,font=2,
        col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clB)])
}

plot_compareClusts_volcano <- function(sCVd,clA,clB,dataType,labType,labNum,labGenes,labAngle) {
  # ^ setup -----
  CGS <- compareClusts_DF(sCVd,clA,clB,dataType)
  CGS <- CGS[!is.na(CGS$FDR),]
  CGS$FDR <- -log10(CGS$FDR)
  temp_exp <- switch(as.character(Param(sCVd,"exponent") == exp(1)),
                     "TRUE"="(natural log scale)",
                     "FALSE"=paste0("(log",Param(sCVd,"exponent")," scale)"))
  if ("flip" %in% labAngle) {
    temp_adjB <- c(-0.1,0.5)
    temp_adjA <- c(1.1,0.5)
  } else {
    temp_adjB <- c(1.1,0.5)
    temp_adjA <- c(-0.1,0.5)
  }
  if (labType == "diff") {
    gnA <- rownames(head(CGS[order(CGS[[dataType]],decreasing=T),],labNum))
    gnB <- rownames(tail(CGS[order(CGS[[dataType]],decreasing=T),],labNum))
  } else if (labType == "de") {
    ts <- order(CGS$FDR,decreasing=T,na.last=T)
    gnA <- rownames(CGS)[ts[CGS[ts,"dir"] == clA][1:labNum]]
    gnB <- rownames(CGS)[ts[CGS[ts,"dir"] == clB][1:labNum]]
  }
  # ^ plot -----
  par(mar=c(3,3,3.5,1),mgp=c(2,1,0))
  plot(x=CGS[[dataType]],y=CGS$FDR,
       xlab=switch(dataType,
                   dDR=paste0("Difference in detection rate (",clA," - ",clB,")"),
                   logGER=paste0("Gene expression ratio (",clA," : ",clB,") ",temp_exp)),
       ylab="-log10 FDR-adjusted p-value",
       main=paste0("Volcano plot of differentially expressed genes (",clA," vs. ",clB,")"),
       pch=20,col=alpha("black",0.3))
  abline(v=0,col="gray50")
  lines(x=c(par("usr")[2],par("usr")[2]),y=c(par("usr")[3],par("usr")[4]),lwd=2,xpd=NA,
        col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clA)])
  lines(x=c(par("usr")[1],par("usr")[1]),y=c(par("usr")[3],par("usr")[4]),xpd=NA,lwd=2,
        col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clB)])
  if (labType == "search") {
    if (length(labGenes) > 0) {
      points(x=CGS[labGenes,dataType],y=CGS[labGenes,"FDR"],
             pch=16,col=alpha("firebrick2",0.8))
      tempLabel <- spreadLabels(x=CGS[labGenes,dataType],y=CGS[labGenes,"FDR"],
                                label=labGenes,str.cex=1.2,str.font=2)
      text(tempLabel,labels=labGenes,col="firebrick2",cex=1.2,font=2)
    }
  } else {
    points(CGS[gnA,dataType],y=CGS[gnA,"FDR"],pch=16,
           col=rainbow2(length(levels(Clusters(sCVd))),.8)[which(levels(Clusters(sCVd)) == clA)])
    points(CGS[gnB,dataType],y=CGS[gnB,"FDR"],pch=16,
           col=rainbow2(length(levels(Clusters(sCVd))),.8)[which(levels(Clusters(sCVd)) == clB)])
    tempLabel <- spreadLabels(CGS[c(gnA,gnB),dataType],CGS[c(gnA,gnB),"FDR"],
                              label=c(gnA,gnB),str.cex=1.2,str.font=2)
    rownames(tempLabel) <- c(gnA,gnB)
    text(tempLabel[gnA,],labels=gnA,cex=1.2,font=2,
         col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clA)])
    text(tempLabel[gnB,],labels=gnB,cex=1.2,font=2,
         col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clB)])
  }
  mtext(paste("Higher in",clA),side=1,line=-1.1,adj=.99,font=2,
        col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clA)])
  mtext(paste("Higher in",clB),side=1,line=-1.1,adj=0.01,font=2,
        col=rainbow2(length(levels(Clusters(sCVd))))[which(levels(Clusters(sCVd)) == clB)])
}


plot_compareClusts <- function(sCVd,clA,clB,dataType,labType,labTypeDiff,labNum,labGenes,labAngle) {
  if (clA %in% levels(Clusters(sCVd)) & 
      clB %in% levels(Clusters(sCVd))) {
    if (dataType %in% c("MGE","MDGE","DR")) {
      plot_compareClusts_MAplot(sCVd,clA,clB,dataType,labType,labNum,labGenes,labAngle)
    } else if (dataType == "GERvDDR") {
      plot_compareClusts_DEscatter(sCVd,clA,clB,dataType,labType,labTypeDiff,labNum,labGenes,labAngle)
    } else if (dataType %in% c("logGER","dDR")) {
      plot_compareClusts_volcano(sCVd,clA,clB,dataType,labType,labNum,labGenes,labAngle)
    }
  } else {
    plot(x=NA,y=NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",xlab=NA,ylab=NA)
    text(.5,.5,paste("Select two clusters to compare in this plot",
                     "using the pulldown menus on the right.",
                     "(Cluster A & Cluster B)",sep="\n"))
  }
}

