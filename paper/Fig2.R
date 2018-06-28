library(pbapply)
library(RColorBrewer)
#### Data setup from PrepareInputs.R ####
mean.logX <- function(data,ex=exponent,pc=pseudocount) { log(mean(ex^data - pc) + 1/ncol(nge),base=ex) }
exponent <- 2  
pseudocount <- 1 

dataRData <- "../scClustViz_files/e13_Cortical_Only.RData" 
convertGeneIDs <- FALSE ##  Set to TRUE if your gene names aren't official gene symbols.

if (exists("dataRDS")) { 
  inD <- readRDS(dataRDS) 
} else if (exists("dataRData")) {
  temp <- load(dataRData)
  inD <- get(temp) ## If you have multiple objects saved in this file, set inD to your data object.
  rm(list=c(temp,"temp"))
} else { warning("Set path to input data as dataRDS or dataRData") }

if (class(inD) == "seurat") {
  require(Seurat)
  inD <- UpdateSeuratObject(inD) ## In case your Seurat object is from an older version of Seurat
  
  nge <- inD@data  
  ##  ^ normalized gene expression matrix (matrix: genes x cells)
  
  if (convertGeneIDs) {
    require(biomaRt)
    e2g <- getBM(attributes=c(geneRowNames,speciesSymbol),
                 mart=mart,filters=geneRowNames,
                 values=rownames(nge))
    e2g <- e2g[e2g[,speciesSymbol] != "",] # removing unmapped gene symbols from conversion table
    print(paste(sum(duplicated(e2g[,geneRowNames])),geneRowNames,"mapped to multiple",speciesSymbol))
    ##  Arbitrarily picking one mapping for the above duplicates, since these generally map to predicted genes anyway.
    e2g <- e2g[!duplicated(e2g[,geneRowNames]),] 
    rownames(e2g) <- e2g[,geneRowNames]
    nge <- nge[e2g[,geneRowNames],] # removing unmapped genes from data
    print(paste(sum(duplicated(e2g[,speciesSymbol])),speciesSymbol,"mapped to multiple",geneRowNames))
    ##  Going to collapse these by summing UMI counts between duplicated rows.
    temp_r <- nge[e2g[,speciesSymbol] %in% e2g[,speciesSymbol][duplicated(e2g[,speciesSymbol])],]
    nge <- nge[!e2g[,speciesSymbol] %in% e2g[,speciesSymbol][duplicated(e2g[,speciesSymbol])],]
    ##  Removed duplicated rows from data, saved as separate object
    rownames(nge) <- e2g[rownames(nge),speciesSymbol] # renamed rows in data as gene symbols
    temp_r <- t(sapply(e2g[,speciesSymbol][duplicated(e2g[,speciesSymbol])],function(X) 
      colSums(temp_r[e2g[,geneRowNames][e2g[,speciesSymbol] == X],])))
    ##  Collapsed by summing each duplicated gene symbol's row
    nge <- rbind(nge,temp_r)  # added those data back to matrix
  }
  
  
  if (!any(grepl("cycle|phase|G2M",colnames(inD@meta.data),ignore.case=T))) {
    data("cc.genes")  
    inD <- CellCycleScoring(inD,g2m.genes=cc.genes$g2m.genes,s.genes=cc.genes$s.genes)
    inD@meta.data$Phase <- factor(inD@meta.data$Phase,levels=c("G1","S","G2M")) # So that the phases are in order.
    rm(cc.genes)
  }
  ##  ^ If necessary, Seurat has a function to predict cell cycle phase from expression of canonical marker genes.
  ##  These are stored as HGNC symbols, but if your data is mouse it will try case-insensitive matches to homologues
  ##  (in which case you will see a warning in AddModuleScore indicating that it attempted to match case).
  
  md <- inD@meta.data[,!grepl("^res",colnames(inD@meta.data))]  
  ##  ^ metadata for cells (dataframe of cells)
  
  if (is.data.frame(inD@meta.data[,grepl("^res",colnames(inD@meta.data))])) {
    cl <- data.frame(lapply(inD@meta.data[,grepl("^res",colnames(inD@meta.data))],as.factor))
  } else {
    cl <- data.frame(inD@meta.data[,grepl("^res",colnames(inD@meta.data))])
    colnames(cl) <- grep("^res",colnames(inD@meta.data),value=T)
  }
  rownames(cl) <- rownames(md) 
  ##  ^ cluster assignments per clustering resolution (dataframe: cells x cluster labels as factors)
  
  if (length(inD@calc.params) == 0) {
    dr_clust <- inD@dr$pca@cell.embeddings
  } else {
    dr_clust <- inD@dr$pca@cell.embeddings[,inD@calc.params$RunTSNE$dims.use]  
  }
  ##  ^ cell embeddings in low-dimensional space used for clustering distances (matrix: cells x dimensions)
  ##  Only including those dimensions used in downstream analysis (ie. those passed to RunTSNE and FindClusters)
  ##  if that information is present (in calc.params).  Else, using all lower dimensions available.
  
  dr_viz <- inD@dr$tsne@cell.embeddings  
  ##  ^ cell embeddings in 2D space for visualization (usually tSNE) (matrix: cells x coordinates)
  
  rm(inD)
} else {
  warning("
          Currently only Seurat objects are supported for auto-loading. 
          You will have to manually copy the data and metadata to the relevant objects. 
          See code above for details."
  )
}

res <- "res.0.8"
CGS <- list()

DR <- pbapply(nge,1,function(X) tapply(X,cl[,res],function(Y) sum(Y>0)/length(Y)))
MDTC <- pbapply(nge,1,function(X) tapply(X,cl[,res],function(Y) {
  temp <- mean.logX(Y[Y>0])
  if (is.na(temp)) { temp <- 0 }
  return(temp)
}))
MTC <- pbapply(nge,1,function(X) tapply(X,cl[,res],mean.logX))
CGS[[res]] <- sapply(levels(cl[,res]),function(X) 
  data.frame(DR=DR[X,],MDTC=MDTC[X,],MTC=MTC[X,]),simplify=F)

combos <- combn(levels(cl[,res]),2)
colnames(combos) <- apply(combos,2,function(X) paste(X,collapse="-"))
deM_dDR <- apply(combos,2,function(i) DR[i[1],] - DR[i[2],])
deM_logFC <- apply(combos,2,function(i) MTC[i[1],] - MTC[i[2],])

pVal_unfiltered <- pbsapply(colnames(combos),function(i)
  apply(nge,1,function(X) 
    wilcox.test(X[cl[,res] == combos[1,i]],
                X[cl[,res] == combos[2,i]])$p.value),simplify=T)
pVal_unfiltered[is.na(pVal_unfiltered)] <- 1


#### Actual experiment ####
par(mar=c(3,3,3,1),mgp=2:0)

plot(as.vector(deM_logFC),-log10(as.vector(pVal_unfiltered)),pch=".",
     ylab="-log10(p-value)",xlab="Log2(Fold Change)")
plot(as.vector(deM_dDR),-log10(as.vector(pVal_unfiltered)),pch=".",
     ylab="-log10(p-value)",xlab=expression(Delta~"Detection Rate"))


TPR <- function(thresh,method) {
  sum(abs(as.vector(method)) > thresh & as.vector(pVal_unfiltered) < 0.01) / sum(as.vector(pVal_unfiltered) < 0.01)
}

FPR <- function(thresh,method) {
  sum(abs(as.vector(method)) > thresh & !(as.vector(pVal_unfiltered) < 0.01)) / sum(!as.vector(pVal_unfiltered) < 0.01)
}

Pr <- function(thresh,method) {
  sum(abs(as.vector(method)) > thresh & as.vector(pVal_unfiltered) < 0.01) / sum(abs(as.vector(method)) > thresh)
}

dDR <- data.frame(TPR=sapply(seq(0,1,.01),TPR,method=deM_dDR),
                  FPR=sapply(seq(0,1,.01),FPR,method=deM_dDR),
                  Pr=sapply(seq(0,1,.01),Pr,method=deM_dDR))
rownames(dDR) <- seq(0,1,.01)
FC <- data.frame(TPR=sapply(seq(0,10,.05),TPR,method=deM_logFC),
                  FPR=sapply(seq(0,10,.05),FPR,method=deM_logFC),
                  Pr=sapply(seq(0,10,.05),Pr,method=deM_logFC))
rownames(FC) <- seq(0,10,.05)

bp <- brewer.pal(3,"Dark2")

plot(TPR~FPR,data=dDR,type="l",col=bp[1],ylab="True Positive Rate",xlab="False Positive Rate",lwd=2)
lines(TPR~FPR,data=FC,type="l",col=bp[2],lwd=2)
legend("top",bty="n",lwd=2,col=bp[1:2],horiz=T,inset=c(0,-.12),xpd=NA,
       legend=c(expression(Delta~"Detection Rate"),"Fold Change"),title="Threshold type")

pdf(file="../Reports/18_F1000Res_scClustViz/Fig2.pdf",width=6,height=6)
par(mar=c(3,3,3,1),mgp=2:0)
plot(Pr~TPR,data=dDR,type="l",col=bp[1],ylab="Precision",xlab="Recall",lwd=2)
lines(Pr~TPR,data=FC,type="l",col=bp[2],lwd=2)
points(Pr~TPR,data=dDR[c("0.1","0.15","0.2"),],pch=19,col=bp[1])
text(dDR[c("0.1","0.15","0.2"),"TPR"],dDR[c("0.1","0.15","0.2"),"Pr"],c("0.1","0.15","0.2"),col=bp[1],pos=4)
legend("top",bty="n",lwd=2,col=bp[1:2],horiz=T,inset=c(0,-.12),xpd=NA,
       legend=c(expression(Delta~"Detection Rate"),"Fold Change"),title="Threshold type")
dev.off()

