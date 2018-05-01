library(pbapply)
library(cluster)

######## User-defined variables ########
exponent <- 2  
##  ^ log base of your normalized input data.  
##  Seurat defaults to natural log (set this to exp(1)), 
##  other methods are generally log2 (set this to 2).
pseudocount <- 1 
##  ^ pseudocount added to all log-normalized values in your input data.  
##  Most methods use a pseudocount of 1 to eliminate log(0) errors.

logFCthresh <- 1 # magnitude of mean log-expression fold change to use as a minimum threshold for DE testing
WRSTalpha <- 0.01 # significance level for DE testing using Wilcoxon rank sum test

dataRDS <- "../scClustViz_files/testData.rds" 
##  ^ path to input data object, saved as RDS (use saveRDS() to generate).
#dataRData <- "testData.RData" 
##  ^ path to input data, saved as RData (use save() to generate )
outputDirectory <- "." 
##  ^ path to output directory (for loading into the R Shiny visualization script)

######## Functions ########
mean.logX <- function(data,ex=exponent,pc=pseudocount) { log(mean(ex^data - pc) + 1/ncol(nge),base=ex) }
##  ^ Adding a pseudocount of 1 to the logMean prior to logFC calculations skews the result quite dramatically,
##  so instead we add a small pseudocount to avoid +/- inf results when means are zero, without the same skewing.
##  Adding a very small (ie 1e-99) number means that means of zero get set to a large negative log-mean, 
##  when it might be more appropriate to have those values fall closer to the smallest non-zero log-mean.
##  By using a pseudocount of 1 / number of samples, we ensure that log(zero) is smaller than any non-zero log-mean,
##  while still being in the same ballpark.
rainbow2 <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 60, c = 100)[1:n]
}

######## Build DE sets for all resolutions ######## 

if (exists("dataRDS")) { 
  inD <- readRDS(dataRDS) 
} else if (exists("dataRData")) {
  temp <- load(dataRData)
  inD <- get(temp)
  rm(list=c(temp,"temp"))
} else { warning("Set path to input data as dataRDS or dataRData") }

if (class(inD) == "seurat") {
  require(Seurat)
  
  nge <- inD@data  
  ##  ^ normalized gene expression matrix (matrix: genes x cells)
  
  if (!any(grepl("cycle|phase|G2M",colnames(inD@meta.data),ignore.case=T))) {
    data("cc.genes")  
    inD <- CellCycleScoring(inD,g2m.genes=cc.genes$g2m.genes,s.genes=cc.genes$s.genes)
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
  
  dr_clust <- inD@dr$pca@cell.embeddings[,inD@calc.params$RunTSNE$dims.use]  
  ##  ^ cell embeddings in low-dimensional space used for clustering distances (matrix: cells x dimensions)
  ##  Only including those dimensions used in downstream analysis (ie. those passed to RunTSNE and FindClusters)
  
  dr_viz <- inD@dr$tsne@cell.embeddings  
  ##  ^ cell embeddings in 2D space for visualization (usually tSNE) (matrix: cells x coordinates)
  
  rm(inD,cc.genes)
} else {
  warning("
Currently only Seurat objects are supported for auto-loading. 
You will have to manually copy the data and metadata to the relevant objects. 
See code above for details."
  )
}

CGS <- deTissue <- deVS <- deMarker <- deNeighb <- list()
for (res in colnames(cl)) {
  #### Precalculate stats for viz tool ####
  print("")
  print("")
  print(paste("Calculating cluster summary statistics for",res))
  print("-- Gene detection rate per cluster --")
  DR <- pbapply(nge,1,function(X) tapply(X,cl[,res],function(Y) sum(Y>0)/length(Y)))
  print("-- Mean detected gene expression per cluster --")
  MDTC <- pbapply(nge,1,function(X) tapply(X,cl[,res],function(Y) {
    temp <- mean.logX(Y[Y>0])
    if (is.na(temp)) { temp <- 0 }
    return(temp)
  }))
  print("-- Mean gene expression per cluster --")
  MTC <- pbapply(nge,1,function(X) tapply(X,cl[,res],mean.logX))
  CGS[[res]] <- sapply(levels(cl[,res]),function(X) 
    data.frame(DR=DR[X,],MDTC=MDTC[X,],MTC=MTC[X,]),simplify=F)
  
  #### deTissue - DE per cluster vs all other data ####
  print("")
  print(paste("Calculating DE vs tissue for",res,"with",length(levels(cl[,res])),"clusters"))
  print("-- LogFC calculations --")
  deT_logFC <- pbsapply(levels(cl[,res]),function(i) 
    MTC[i,] - apply(nge[,cl[,res] != i],1,mean.logX))
  deT_genesUsed <- apply(deT_logFC,2,function(X) which(X > logFCthresh))  
  if (any(sapply(deT_genesUsed,length) < 1)) {
    stop(paste0("logFCthresh should be set to less than ",
                min(apply(deT_logFC,2,function(X) max(abs(X)))),
                ", the largest magnitude logFC between cluster ",
                names(which.min(apply(deT_logFC,2,function(X) max(abs(X))))),
                " and the remaining data."))
  }
  print("-- Wilcoxon rank sum calculations --")
  deT_pVal <- pbsapply(levels(cl[,res]),function(i)
    apply(nge[deT_genesUsed[[i]],],1,function(X) 
      wilcox.test(X[cl[,res] == i],X[cl[,res] != i])$p.value),simplify=F)
  deTissue[[res]] <- sapply(levels(cl[,res]),function(i) 
    data.frame(logFC=deT_logFC[deT_genesUsed[[i]],i],
               pVal=deT_pVal[[i]])[order(deT_pVal[[i]]),],simplify=F)
  tempQval <- tapply(p.adjust(do.call(rbind,deTissue[[res]])$pVal,"fdr"),
                     rep(names(sapply(deTissue[[res]],nrow)),sapply(deTissue[[res]],nrow)),c)
  for (i in names(deTissue[[res]])) { deTissue[[res]][[i]]$qVal <- tempQval[[i]] }
  
  #### deMarker - DE per cluster vs each other cluster #### 
  combos <- combn(levels(cl[,res]),2)
  colnames(combos) <- apply(combos,2,function(X) paste(X,collapse="-"))
  print("")
  print(paste("Calculating marker DE for",res,"with",ncol(combos),"combinations of clusters"))
  deM_dDR <- apply(combos,2,function(i) DR[i[1],] - DR[i[2],])
  deM_logFC <- apply(combos,2,function(i) MTC[i[1],] - MTC[i[2],])
  deM_genesUsed <- apply(deM_logFC,2,function(X) which(abs(X) > logFCthresh))  
  if (any(sapply(deM_genesUsed,length) < 1)) {
    stop(paste0("logFCthresh should be set to less than ",
                min(apply(deM_logFC,2,function(X) max(abs(X)))),
                ", the largest magnitude logFC between clusters ",
                names(which.min(apply(deM_logFC,2,function(X) max(abs(X))))),"."))
  }
  deM_pVal <- pbsapply(colnames(combos),function(i)
    apply(nge[deM_genesUsed[[i]],],1,function(X) 
      wilcox.test(X[cl[,res] == combos[1,i]],
                  X[cl[,res] == combos[2,i]])$p.value),simplify=F)
  temp_deVS <- sapply(colnames(combos),function(i) 
    data.frame(dDR=deM_dDR[deM_genesUsed[[i]],i],logFC=deM_logFC[deM_genesUsed[[i]],i],
               pVal=deM_pVal[[i]])[order(deM_pVal[[i]]),],simplify=F)
  tempQval <- tapply(p.adjust(do.call(rbind,temp_deVS)$pVal,"fdr"),
                     rep(names(sapply(temp_deVS,nrow)),sapply(temp_deVS,nrow)),c)
  for (i in names(temp_deVS)) { temp_deVS[[i]]$qVal <- tempQval[[i]] }
  
  deVS[[res]] <- sapply(levels(cl[,res]),function(i) {
    combos <- strsplit(names(temp_deVS),"-")
    temp <- list()
    for (X in seq_along(combos)) {
      if (! i %in% combos[[X]]) {
        next
      } else if (which(combos[[X]] == i) == 1) {
        temp[[combos[[X]][2]]] <- temp_deVS[[X]][temp_deVS[[X]]$logFC > 0 & temp_deVS[[X]]$qVal < WRSTalpha,]
      } else if (which(combos[[X]] == i) == 2) {
        temp[[combos[[X]][1]]] <- temp_deVS[[X]][temp_deVS[[X]]$logFC < 0 & temp_deVS[[X]]$qVal < WRSTalpha,]
        temp[[combos[[X]][1]]]$dDR <- temp[[combos[[X]][1]]]$dDR * -1
        temp[[combos[[X]][1]]]$logFC <- temp[[combos[[X]][1]]]$logFC * -1
      }
    }
    return(temp)
  },simplify=F)
  
  deMarker[[res]] <- sapply(deVS[[res]],function(X) {
    markerGenes <- Reduce(intersect,lapply(X,rownames))
    temp <- sapply(X,function(Y) Y[markerGenes,c("dDR","logFC","qVal")],simplify=F)
    names(temp) <- paste("vs",names(temp),sep=".")
    return(do.call(cbind,temp))
  },simplify=F)
  
  ### deNeighb - DE between closest neighbouring clusters ####
  nb <- apply(dist(apply(dr_viz,2,
                         function(X) tapply(X,cl[,res],mean)),diag=T,upper=T),2,
              function(Z) names(which.min(Z[Z > 0])))
  
  deNeighb[[res]] <- mapply(function(NB,VS) VS[[NB]][,c("dDR","logFC","qVal")],NB=nb,VS=deVS[[res]],SIMPLIFY=F)
  for (i in names(deNeighb[[res]])) {
    colnames(deNeighb[[res]][[i]]) <- paste("vs",nb[i],colnames(deNeighb[[res]][[i]]),sep=".")
  }
}

#### Save outputs for visualization ####
save(nge,md,cl,dr_clust,dr_viz,
     CGS,deTissue,deVS,deMarker,deNeighb,
     file=paste0(sub("\\..+$","",get(grep("^dataRD",ls(),value=T))),"_forViz.RData"))
##  ^ Saved objects for use in visualization script (RunVizScript.R).
