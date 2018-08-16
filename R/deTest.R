#' Cluster-wise differential expression testing
#'
#' Performs differential expression testing between clusters for all cluster solutions in 
#' order to assess the biological relevance of each cluster solution. Differential 
#' expression testing is done using the Wilcoxon rank-sum test implemented in the base R 
#' \code{stats} package. For details about what is being compared in the tests, see the 
#' "Value" section.
#'
#' @param il The list outputted by one of the importData functions (either 
#'   \code{\link{readFromSeurat}} or \code{\link{readFromManual}}).
#'   
#' @param testAll Logical value indicating whether to test all cluster solutions 
#'   (\code{TRUE}) or stop testing once a cluster solution has been found where there is 
#'   no differentially expressed genes found between at least one pair of nearest 
#'   neighbouring clusters (\code{FALSE}). \emph{If set to (\code{FALSE}), only the 
#'   cluster solutions tested will appear in the scClustViz shiny app.}
#' 
#' @param exponent The log base of your normalized input data. Seurat normalization uses 
#'   the natural log (set this to exp(1)), while other normalization methods generally use 
#'   log2 (set this to 2).
#' 
#' @param pseudocount The pseudocount added to all log-normalized values in your input 
#'   data. Most methods use a pseudocount of 1 to eliminate log(0) errors.
#' 
#' @param FDRthresh The false discovery rate to use as a threshold for determining statistical 
#'   significance of differential expression calculated by the Wilcoxon rank-sum test.
#' 
#' @param threshType Filtering genes for use in differential expression testing can be 
#'   done multiple ways. We use an expression ratio filter for comparing each cluster to 
#'   the rest of the tissue as a whole, but find that difference in detection rates works 
#'   better when comparing clusters to each other. You can set threshType to 
#'   \code{"logGER"} to use a gene expression ratio for all gene filtering, or leave it as 
#'   default (\code{"dDR"}) to use difference in detection rate as the thresholding method 
#'   when comparing clusters to each other.
#' 
#' @param dDRthresh Magnitude of detection rate difference of a gene between clusters to 
#'   use as filter for determining which genes to test for differential expression between 
#'   clusters.
#' 
#' @param logGERthresh Magnitude of gene expression ratio for a gene between clusters to 
#'   use as filter for determining which genes to test for differential expression between 
#'   clusters.
#' 
#' @return The function returns a list containing the results of differential expression 
#'   testing for all sets of cluster solutions. \emph{Saving both the input (the object passed 
#'   to the \code{il} argument) and the output of this function to an RData file is all 
#'   the preparation necessary for running the scClustViz Shiny app itself.}
#'   The output list of this function contains the following elements:
#'   \describe{
#'     \item{CGS}{} 
#'     \item{deTissue}{} 
#'     \item{deVS}{}
#'     \item{deMarker}{}
#'     \item{deDist}{} 
#'     \item{deNeighb}{} 
#'   }
#' 
#' @examples 
#' \dontrun{
#'  data_for_scClustViz <- readFromSeurat(your_seurat_object,
#'                                        convertGeneIDs=F)
#'  rm(your_seurat_object) 
#'  # All the data scClustViz needs is in 'data_for_scClustViz'.
#'  
#'  DE_for_scClustViz <- clusterWiseDEtest(data_for_scClustViz)
#'  
#'  save(data_for_scClustViz,DE_for_scClustViz,
#'       file="for_scClustViz.RData")
#'  # Save these objects so you'll never have to run this slow function again!
#'  
#'  runShiny(filePath="for_scClustViz.RData")
#' }
#' 
#' @seealso \code{\link{readFromSeurat}} or \code{\link{readFromManual}} for reading in 
#'   data to generate the input object for this function, and \code{\link{runShiny}} to 
#'   use the interactive Shiny GUI to view the results of this testing.
#' 
#' @export




clusterWiseDEtest <- function(il,testAll=TRUE,
                              exponent=2,pseudocount=1,FDRthresh=0.01,
                              threshType="dDR",dDRthresh=0.15,logGERthresh=1) {
  temp_warn <- options("warn")
  options(warn=-1)
  
  out <- list(CGS=list(),deTissue=list(),deVS=list(),
              deMarker=list(),deDist=list(),deNeighb=list())
  # This loop iterates through every cluster solution, and does DE testing between clusters
  # to generate the DE metrics for assessing your clusters.  This takes some time.
  for (res in colnames(il[["cl"]])) {
    #### Precalculate stats for viz tool ####
    print("")
    print("")
    print(paste("Calculating cluster gene summary statistics for",res))
    print("-- Gene detection rate per cluster --")
    DR <- pbapply::pbapply(il[["nge"]],1,function(X) tapply(X,il[["cl"]][,res],function(Y) sum(Y>0)/length(Y)))
    print("-- Mean detected gene expression per cluster --")
    MDTC <- pbapply::pbapply(il[["nge"]],1,function(X) tapply(X,il[["cl"]][,res],function(Y) {
      temp <- meanLogX(Y[Y>0],ncell=ncol(il[["nge"]]),ex=exponent,pc=pseudocount)
      if (is.na(temp)) { temp <- 0 }
      return(temp)
    }))
    print("-- Mean gene expression per cluster --")
    MTC <- pbapply::pbapply(il[["nge"]],1,function(X) 
      tapply(X,il[["cl"]][,res],function(Y) 
        meanLogX(Y,ncell=ncol(il[["nge"]]),ex=exponent,pc=pseudocount)))
    out[["CGS"]][[res]] <- sapply(levels(il[["cl"]][,res]),function(X) 
      data.frame(DR=DR[X,],MDTC=MDTC[X,],MTC=MTC[X,]),simplify=F)
    
    #### deTissue - DE per cluster vs all other data ####
    print("")
    print(paste("Calculating DE vs tissue for",res,"with",length(levels(il[["cl"]][,res])),"clusters"))
    print("-- logGER calculations --")
    deT_logGER <- pbapply::pbsapply(levels(il[["cl"]][,res]),function(i) 
      MTC[i,] - apply(il[["nge"]][,il[["cl"]][,res] != i],1,function(Y) 
        meanLogX(Y,ncell=ncol(il[["nge"]]),ex=exponent,pc=pseudocount)))
    deT_genesUsed <- apply(deT_logGER,2,function(X) which(X > logGERthresh))  
    if (any(sapply(deT_genesUsed,length) < 1)) {
      stop(paste0("logGERthresh should be set to less than ",
                  min(apply(deT_logGER,2,function(X) max(abs(X)))),
                  ", the largest magnitude logGER between cluster ",
                  names(which.min(apply(deT_logGER,2,function(X) max(abs(X))))),
                  " and the remaining data."))
    }
    print("-- Wilcoxon rank sum calculations --")
    deT_pVal <- pbapply::pbsapply(levels(il[["cl"]][,res]),function(i)
      apply(il[["nge"]][deT_genesUsed[[i]],],1,function(X) 
        wilcox.test(X[il[["cl"]][,res] == i],X[il[["cl"]][,res] != i])$p.value),simplify=F)
    out[["deTissue"]][[res]] <- sapply(levels(il[["cl"]][,res]),function(i) 
      data.frame(logGER=deT_logGER[deT_genesUsed[[i]],i],
                 pVal=deT_pVal[[i]])[order(deT_pVal[[i]]),],simplify=F)
    tempQval <- tapply(p.adjust(do.call(rbind,out[["deTissue"]][[res]])$pVal,"fdr"),
                       rep(names(sapply(out[["deTissue"]][[res]],nrow)),sapply(out[["deTissue"]][[res]],nrow)),c)
    for (i in names(out[["deTissue"]][[res]])) { 
      out[["deTissue"]][[res]][[i]] <- out[["deTissue"]][[res]][[i]][tempQval[[i]] <= FDRthresh,]
      out[["deTissue"]][[res]][[i]]$qVal <- tempQval[[i]][tempQval[[i]] <= FDRthresh] 
    }
    
    #### deMarker - DE per cluster vs each other cluster #### 
    combos <- combn(levels(il[["cl"]][,res]),2)
    colnames(combos) <- apply(combos,2,function(X) paste(X,collapse="-"))
    print("")
    print(paste("Calculating marker DE for",res,"with",ncol(combos),"combinations of clusters"))
    deM_dDR <- apply(combos,2,function(i) DR[i[1],] - DR[i[2],])
    deM_logGER <- apply(combos,2,function(i) MTC[i[1],] - MTC[i[2],])
    deM_genesUsed <- switch(threshType,
                            dDR=apply(deM_dDR,2,function(X) which(abs(X) > dDRthresh)),
                            logGER=apply(deM_logGER,2,function(X) which(abs(X) > logGERthresh))) 
    if (any(sapply(deM_genesUsed,length) < 1)) {
      stop("Gene filtering threshold is set too high.")
    }
    deM_pVal <- pbapply::pbsapply(colnames(combos),function(i)
      apply(il[["nge"]][deM_genesUsed[[i]],],1,function(X) 
        wilcox.test(X[il[["cl"]][,res] == combos[1,i]],
                    X[il[["cl"]][,res] == combos[2,i]])$p.value),simplify=F)
    temp_deVS <- sapply(colnames(combos),function(i) 
      data.frame(dDR=deM_dDR[deM_genesUsed[[i]],i],logGER=deM_logGER[deM_genesUsed[[i]],i],
                 pVal=deM_pVal[[i]])[order(deM_pVal[[i]]),],simplify=F)
    tempQval <- tapply(p.adjust(do.call(rbind,temp_deVS)$pVal,"fdr"),
                       rep(names(sapply(temp_deVS,nrow)),sapply(temp_deVS,nrow)),c)
    for (i in names(temp_deVS)) { temp_deVS[[i]]$qVal <- tempQval[[i]] }
    
    out[["deVS"]][[res]] <- sapply(levels(il[["cl"]][,res]),function(i) {
      combos <- strsplit(names(temp_deVS),"-")
      temp <- list()
      for (X in seq_along(combos)) {
        if (! i %in% combos[[X]]) {
          next
        } else if (which(combos[[X]] == i) == 1) {
          temp[[combos[[X]][2]]] <- temp_deVS[[X]][temp_deVS[[X]][,threshType] > 0 & 
                                                     temp_deVS[[X]]$qVal <= FDRthresh,]
        } else if (which(combos[[X]] == i) == 2) {
          temp[[combos[[X]][1]]] <- temp_deVS[[X]][temp_deVS[[X]][,threshType] < 0 &
                                                     temp_deVS[[X]]$qVal <= FDRthresh,]
          temp[[combos[[X]][1]]]$dDR <- temp[[combos[[X]][1]]]$dDR * -1
          temp[[combos[[X]][1]]]$logGER <- temp[[combos[[X]][1]]]$logGER * -1
        }
      }
      return(temp)
    },simplify=F)
    
    out[["deMarker"]][[res]] <- sapply(out[["deVS"]][[res]],function(X) {
      markerGenes <- Reduce(intersect,lapply(X,rownames))
      temp <- sapply(X,function(Y) Y[markerGenes,c("dDR","logGER","qVal")],simplify=F)
      names(temp) <- paste("vs",names(temp),sep=".")
      return(do.call(cbind,temp))
    },simplify=F)
    
    ### deNeighb - DE between nearest neighbouring clusters ####
    out[["deDist"]][[res]] <- sapply(names(out[["deVS"]][[res]]),function(X) sapply(names(out[["deVS"]][[res]]),function(Y) 
      if (X == Y) { return(NA) } else { min(nrow(out[["deVS"]][[res]][[X]][[Y]]),nrow(out[["deVS"]][[res]][[Y]][[X]])) }))
    nb <- colnames(out[["deDist"]][[res]])[apply(out[["deDist"]][[res]],1,which.min)]
    names(nb) <- colnames(out[["deDist"]][[res]])
    ##  Nearest neighbour determined by number of DE genes between clusters.
    
    out[["deNeighb"]][[res]] <- mapply(function(NB,VS) VS[[NB]][,c("dDR","logGER","qVal")],
                                       NB=nb,VS=out[["deVS"]][[res]],SIMPLIFY=F)
    for (i in names(out[["deNeighb"]][[res]])) {
      colnames(out[["deNeighb"]][[res]][[i]]) <- paste("vs",nb[i],colnames(out[["deNeighb"]][[res]][[i]]),sep=".")
    }
    
    if (!testAll) {
      if (min(sapply(out[["deNeighb"]][[res]],nrow)) < 1) { break }
    }
  }
  options(warn=temp_warn$warn)
  return(out)
}
