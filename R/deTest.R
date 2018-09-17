#' Cluster-wise differential expression testing
#'
#' Performs differential expression testing between clusters for all cluster
#' solutions in order to assess the biological relevance of each cluster
#' solution. Differential expression testing is done using the Wilcoxon rank-sum
#' test implemented in the base R \code{stats} package. For details about what
#' is being compared in the tests, see the "Value" section. This is a wrapper 
#' function for running \code{\link{calcAllDE}} over each cluster resolution in 
#' the input.
#'
#' @param il The list outputted by one of the importData functions (either
#'   \code{\link{readFromSeurat}} or \code{\link{readFromManual}}).
#'
#' @param testAll Default = TRUE. Logical value indicating whether to test all
#'   cluster solutions (\code{TRUE}) or stop testing once a cluster solution has
#'   been found where there is no differentially expressed genes found between
#'   at least one pair of nearest neighbouring clusters (\code{FALSE}). \emph{If
#'   set to (\code{FALSE}), only the cluster solutions tested will appear in the
#'   scClustViz shiny app.}
#'
#' @param exponent Default = 2. The log base of your normalized input data.
#'   Seurat normalization uses the natural log (set this to exp(1)), while other
#'   normalization methods generally use log2 (set this to 2).
#'
#' @param pseudocount Default = 1. The pseudocount added to all log-normalized
#'   values in your input data. Most methods use a pseudocount of 1 to eliminate
#'   log(0) errors.
#'
#' @param FDRthresh Default = 0.01. The false discovery rate to use as a
#'   threshold for determining statistical significance of differential
#'   expression calculated by the Wilcoxon rank-sum test.
#'
#' @param threshType Default = "dDR". Filtering genes for use in differential
#'   expression testing can be done multiple ways. We use an expression ratio
#'   filter for comparing each cluster to the rest of the tissue as a whole, but
#'   find that difference in detection rates works better when comparing
#'   clusters to each other. You can set threshType to \code{"logGER"} to use a
#'   gene expression ratio for all gene filtering, or leave it as default
#'   (\code{"dDR"}) to use difference in detection rate as the thresholding
#'   method when comparing clusters to each other.
#'
#' @param dDRthresh Default = 0.15. Magnitude of detection rate difference of a
#'   gene between clusters to use as filter for determining which genes to test
#'   for differential expression between clusters.
#'
#' @param logGERthresh Default = 1. Magnitude of gene expression ratio for a
#'   gene between clusters to use as filter for determining which genes to test
#'   for differential expression between clusters.
#'
#' @return The function returns a list containing the results of differential
#'   expression testing for all sets of cluster solutions. \emph{Saving both the
#'   input (the object passed to the \code{il} argument) and the output of this
#'   function to an RData file is all the preparation necessary for running the
#'   scClustViz Shiny app itself.} The output list of this function contains the
#'   following elements: 
#'   \describe{ 
#'     \item{CGS}{A nested list of dataframes. Each list element is named for 
#'       a column in \code{il$cl} (a cluster resolution). That list element 
#'       contains a named list of clusters at that resolution. Each of those 
#'       list elements contains a dataframe of three variables, where each 
#'       sample is a gene. \code{DR} is the proportion of cells in the cluster 
#'       in which that gene was detected. \code{MDTC} is mean normalized gene 
#'       expression for that gene in only the cells in which it was detected 
#'       (see \link{meanLogX} for mean calculation). \code{MTC} is the mean 
#'       normalized gene expression for that gene in all cells of the cluster 
#'       (see \link{meanLogX} for mean calculation).} 
#'     \item{deTissue}{Differential testing results from Wilcoxon rank sum tests 
#'       comparing a gene in each cluster to the rest of the cells as a whole in 
#'       a one vs all comparison. The results are stored as a nested list of 
#'       dataframes. Each list element is named for a column in \code{il$cl} (a 
#'       cluster resolution). That list element contains a named list of 
#'       clusters at that resolution. Each of those list elements contains a 
#'       dataframe of three variables, where each sample is a gene. 
#'       \code{logGER} is the log gene expression ratio calculated by 
#'       subtracting the mean expression of the gene (see \link{meanLogX} for 
#'       mean calculation) in all other cells from the mean expression of the 
#'       gene in this cluster. \code{pVal} is the p-value of the Wilcoxon rank 
#'       sum test. \code{qVal} is the false discovery rate-corrected p-value of 
#'       the test.} 
#'     \item{deVS}{Differential testing results from Wilcoxon rank sum tests 
#'       comparing a gene in each cluster to that gene in every other cluster in 
#'       a series of tests. The results are stored as a nested list of 
#'       dataframes. Each list element is named for a column in \code{il$cl} (a 
#'       cluster resolution). That list element contains a named list of 
#'       clusters at that resolution (cluster A). Each of those lists contains a 
#'       named list of all the other clusters at that resolution (cluster B). 
#'       Each of those list elements contains a dataframe of four variables, 
#'       where each sample is a gene. \code{dDR} is the difference in detection 
#'       rate of that gene between the two clusters (DR[A] - DR[B]). 
#'       \code{logGER} is the log gene expression ratio calculated by taking the 
#'       difference in mean expression of the gene (see \link{meanLogX} for 
#'       mean calculation) between the two clusters (MTC[A] - MTC[B]). 
#'       \code{pVal} is the p-value of the Wilcoxon rank sum test. \code{qVal} 
#'       is the false discovery rate-corrected p-value of the test.}
#'     \item{deMarker}{Differential testing results from Wilcoxon rank sum tests 
#'       comparing a gene in each cluster to that gene in every other cluster in 
#'       a series of tests, and filtering for only those genes that show 
#'       significant positive differential expression versus all other clusters. 
#'       The results are stored as a nested list of dataframes. Each list 
#'       element is named for a column in \code{il$cl} (a cluster resolution). 
#'       That list element contains a named list of clusters at that resolution 
#'       (cluster A). Each of those list elements contains a dataframe where 
#'       variables represent comparisons to all the other clusters and each 
#'       sample is a gene. For each other cluster (cluster B), there are three 
#'       variables, named as follows: \code{vs.B.dDR} is the difference in 
#'       detection rate of that gene between the two clusters (DR[A] - DR[B]). 
#'       \code{vs.B.logGER} is the log gene expression ratio calculated by 
#'       taking the difference in mean expression of the gene (see 
#'       \link{meanLogX} for mean calculation) between the two clusters (MTC[A] 
#'       - MTC[B]). \code{vs.B.qVal} is the false discovery rate-corrected 
#'       p-value of the Wilcoxon rank sum test.} 
#'     \item{deDist}{A named list of distances between clusters for each cluster 
#'       resolution. Distances are calculated as number of differentially 
#'       expressed genes between clusters.} 
#'     \item{deNeighb}{Differential testing results from Wilcoxon rank sum tests 
#'       comparing a gene in each cluster to that gene in its nearest 
#'       neighbouring cluster (calculated by number of differentially expressed 
#'       genes), and filtering for only those genes that show significant 
#'       positive differential expression versus all other clusters. The results 
#'       are stored as a nested list of dataframes. Each list element is named 
#'       for a column in \code{il$cl} (a cluster resolution). That list element 
#'       contains a named list of clusters at that resolution (cluster A). Each 
#'       of those list elements contains a dataframe where variables represent 
#'       the comparison to its nearest neighbouring cluster (cluster B) and each 
#'       sample is a gene. There are three variables, named as follows: 
#'       \code{vs.B.dDR} is the difference in detection rate of that gene 
#'       between the two clusters (DR[A] - DR[B]). \code{vs.B.logGER} is the log 
#'       gene expression ratio calculated by taking the difference in mean 
#'       expression of the gene (see \link{meanLogX} for mean calculation) 
#'       between the two clusters (MTC[A] - MTC[B]). \code{vs.B.qVal} is the 
#'       false discovery rate-corrected p-value of the Wilcoxon rank sum test.}
#'     \item{params}{A list of the parameters from the argument list of this 
#'       function used to do the analysis, saved so that the same parameters are 
#'       used in the Shiny app.} 
#'   }
#'
#' @examples
#' \dontrun{
#'  data_for_scClustViz <- readFromSeurat(your_seurat_object)
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
#' @seealso \code{\link{readFromSeurat}} or \code{\link{readFromManual}} for
#'   reading in data to generate the input object for this function, and
#'   \code{\link{runShiny}} to use the interactive Shiny GUI to view the results
#'   of this testing.
#'
#' @export

clusterWiseDEtest <- function(il,testAll=TRUE,
                              exponent=2,pseudocount=1,FDRthresh=0.01,
                              threshType="dDR",dDRthresh=0.15,logGERthresh=1) {
  
  out <- list(CGS=list(),deTissue=list(),deVS=list(),
              deMarker=list(),deDist=list(),deNeighb=list(),
              params=list(exponent=exponent,
                          pseudocount=pseudocount,
                          FDRthresh=FDRthresh,
                          threshType=threshType,
                          dDRthresh=dDRthresh,
                          logGERthresh=logGERthresh))
  # This loop iterates through every cluster solution, and does DE testing between clusters
  # to generate the DE metrics for assessing your clusters.  This takes some time.
  for (res in colnames(il[["cl"]])) {
    print("")
    print("")
    print(paste("Calculating all DE stats for ",res))
    tempOut <- calcAllDE(nge=il[["nge"]],
                         cl=il[["cl"]][[res]],
                         params=out$params)
    out$CGS[[res]] <- tempOut$CGS
    out$deTissue[[res]] <- tempOut$deTissue
    out$deVS[[res]] <- tempOut$deVS
    out$deMarker[[res]] <- tempOut$deMarker
    out$deDist[[res]] <- tempOut$deDist
    out$deNeighb[[res]] <- tempOut$deNeighb
    if (!testAll) {
      if (min(sapply(tempOut$deNeighb,nrow)) < 1) { break }
    }
  }
  return(out)
}


#' Calculate gene stats and differential expression for a single cluster resolution
#'
#' Performs differential expression testing between clusters in order to assess 
#' the biological relevance of each cluster solution. Differential expression 
#' testing is done using the Wilcoxon rank-sum test implemented in the base R 
#' \code{stats} package. For details about what is being compared in the tests, 
#' see the "Value" section. You probably don't need to use this unless you're 
#' incorporating \code{scClustViz} into your analysis pipeline and testing DE 
#' after each clustering run. Otherwise see \code{\link{clusterWiseDEtest}}.
#'
#' @param nge The log-normalized gene expression matrix.
#'
#' @param cl The factor with cluster assignments per cell (column of nge).
#'
#' @param params A list with the following parameters:
#'   \describe{
#'     \item{exponent}{The log base of your normalized input data.
#'       Seurat normalization uses the natural log (set this to exp(1)), while other
#'       normalization methods generally use log2 (set this to 2).}
#'     \item{pseudocount}{The pseudocount added to all log-normalized
#'       values in your input data. Most methods use a pseudocount of 1 to eliminate
#'       log(0) errors.}
#'     \item{FDRthresh}{The false discovery rate to use as a
#'       threshold for determining statistical significance of differential
#'       expression calculated by the Wilcoxon rank-sum test.}
#'     \item{threshType}{Filtering genes for use in differential
#'       expression testing can be done multiple ways. We use an expression ratio
#'       filter for comparing each cluster to the rest of the tissue as a whole, but
#'       find that difference in detection rates works better when comparing
#'       clusters to each other. You can set threshType to \code{"logGER"} to use a
#'       gene expression ratio for all gene filtering, or leave it as default
#'       (\code{"dDR"}) to use difference in detection rate as the thresholding
#'       method when comparing clusters to each other.}
#'     \item{dDRthresh}{Magnitude of detection rate difference of a
#'       gene between clusters to use as filter for determining which genes to test
#'       for differential expression between clusters.}
#'     \item{logGERthresh}{Magnitude of gene expression ratio for a
#'        gene between clusters to use as filter for determining which genes to test
#'        for differential expression between clusters.}
#'   } 
#'
#' @return The function returns a list containing the results of differential
#'   expression testing for a cluster solution. The output list of this function 
#'   contains the following elements: 
#'   \describe{ 
#'     \item{CGS}{A nested list of dataframes. Each list element 
#'       contains a named list of clusters. Each of those 
#'       list elements contains a dataframe of three variables, where each 
#'       sample is a gene. \code{DR} is the proportion of cells in the cluster 
#'       in which that gene was detected. \code{MDTC} is mean normalized gene 
#'       expression for that gene in only the cells in which it was detected 
#'       (see \link{meanLogX} for mean calculation). \code{MTC} is the mean 
#'       normalized gene expression for that gene in all cells of the cluster 
#'       (see \link{meanLogX} for mean calculation).} 
#'     \item{deTissue}{Differential testing results from Wilcoxon rank sum tests 
#'       comparing a gene in each cluster to the rest of the cells as a whole in 
#'       a one vs all comparison. The results are stored as a nested list of 
#'       dataframes. Each list element contains a named list of 
#'       clusters. Each of those list elements contains a 
#'       dataframe of three variables, where each sample is a gene. 
#'       \code{logGER} is the log gene expression ratio calculated by 
#'       subtracting the mean expression of the gene (see \link{meanLogX} for 
#'       mean calculation) in all other cells from the mean expression of the 
#'       gene in this cluster. \code{pVal} is the p-value of the Wilcoxon rank 
#'       sum test. \code{qVal} is the false discovery rate-corrected p-value of 
#'       the test.} 
#'     \item{deVS}{Differential testing results from Wilcoxon rank sum tests 
#'       comparing a gene in each cluster to that gene in every other cluster in 
#'       a series of tests. The results are stored as a nested list of 
#'       dataframes. Each list element contains a named list of 
#'       clusters (cluster A). Each of those lists contains a 
#'       named list of all the other clusters at that resolution (cluster B). 
#'       Each of those list elements contains a dataframe of four variables, 
#'       where each sample is a gene. \code{dDR} is the difference in detection 
#'       rate of that gene between the two clusters (DR[A] - DR[B]). 
#'       \code{logGER} is the log gene expression ratio calculated by taking the 
#'       difference in mean expression of the gene (see \link{meanLogX} for 
#'       mean calculation) between the two clusters (MTC[A] - MTC[B]). 
#'       \code{pVal} is the p-value of the Wilcoxon rank sum test. \code{qVal} 
#'       is the false discovery rate-corrected p-value of the test.}
#'     \item{deMarker}{Differential testing results from Wilcoxon rank sum tests 
#'       comparing a gene in each cluster to that gene in every other cluster in 
#'       a series of tests, and filtering for only those genes that show 
#'       significant positive differential expression versus all other clusters. 
#'       The results are stored as a nested list of dataframes. Each list element 
#'       contains a named list of clusters (cluster A). Each of those list 
#'       elements contains a dataframe where 
#'       variables represent comparisons to all the other clusters and each 
#'       sample is a gene. For each other cluster (cluster B), there are three 
#'       variables, named as follows: \code{vs.B.dDR} is the difference in 
#'       detection rate of that gene between the two clusters (DR[A] - DR[B]). 
#'       \code{vs.B.logGER} is the log gene expression ratio calculated by 
#'       taking the difference in mean expression of the gene (see 
#'       \link{meanLogX} for mean calculation) between the two clusters (MTC[A] 
#'       - MTC[B]). \code{vs.B.qVal} is the false discovery rate-corrected 
#'       p-value of the Wilcoxon rank sum test.} 
#'     \item{deDist}{A matrix of distances between clusters.
#'       Distances are calculated as number of differentially 
#'       expressed genes between clusters.} 
#'     \item{deNeighb}{Differential testing results from Wilcoxon rank sum tests 
#'       comparing a gene in each cluster to that gene in its nearest 
#'       neighbouring cluster (calculated by number of differentially expressed 
#'       genes), and filtering for only those genes that show significant 
#'       positive differential expression versus all other clusters. The results 
#'       are stored as a nested list of dataframes. Each list element 
#'       contains a named list of clusters (cluster A). Each 
#'       of those list elements contains a dataframe where variables represent 
#'       the comparison to its nearest neighbouring cluster (cluster B) and each 
#'       sample is a gene. There are three variables, named as follows: 
#'       \code{vs.B.dDR} is the difference in detection rate of that gene 
#'       between the two clusters (DR[A] - DR[B]). \code{vs.B.logGER} is the log 
#'       gene expression ratio calculated by taking the difference in mean 
#'       expression of the gene (see \link{meanLogX} for mean calculation) 
#'       between the two clusters (MTC[A] - MTC[B]). \code{vs.B.qVal} is the 
#'       false discovery rate-corrected p-value of the Wilcoxon rank sum test.}
#'   }
#'
#' @examples
#' \dontrun{
#'  data_for_scClustViz <- readFromSeurat(your_seurat_object)
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
#' @seealso \code{\link{readFromSeurat}} or \code{\link{readFromManual}} for
#'   reading in data to generate the input object for this function, and
#'   \code{\link{runShiny}} to use the interactive Shiny GUI to view the results
#'   of this testing.
#'
#' @export

calcAllDE <- function(nge,cl,params) {
  out <- list()
  print("")
  print("Calculating cluster gene summary statistics")
  out[["CGS"]] <- calcCGS(nge=nge,
                          cl=cl,
                          exponent=params$exponent,
                          pseudocount=params$pseudocount)
  print("")
  print(paste("Calculating DE vs tissue with",length(levels(cl)),"clusters"))
  print("-- logGER calculations --")
  deTes <- calcESvsRest(nge=nge,
                        cl=cl,
                        CGS=out[["CGS"]],
                        exponent=params$exponent,
                        pseudocount=params$pseudocount,
                        logGERthresh=params$logGERthresh)
  if (any(sapply(deTes,length) < 1)) {
    stop("Gene filtering threshold is set too high.")
  }
  print("-- Wilcoxon rank sum calculations --")
  out[["deTissue"]] <- calcDEvsRest(nge,
                                    cl=cl,
                                    deTes=deTes,
                                    exponent=params$exponent,
                                    pseudocount=params$pseudocount,
                                    FDRthresh=params$FDRthresh)
  print("")
  print(paste("Calculating marker DE with",ncol(combn(levels(cl),2)),"combinations of clusters"))
  deMes <- calcESvsCombn(cl=cl,
                         CGS=out[["CGS"]],
                         threshType=params$threshType,
                         logGERthresh=params$logGERthresh,
                         dDRthresh=params$dDRthresh)
  if (any(sapply(deMes,nrow) < 1)) {
    stop("Gene filtering threshold is set too high.")
  }
  out[["deVS"]] <- calcDEvsCombn(nge=nge,
                                 cl=cl,
                                 deMes=deMes,
                                 threshType=params$threshType,
                                 FDRthresh=params$FDRthresh)
  out[["deMarker"]] <- calcMarker(deVS=out[["deVS"]],
                                  FDRthresh=params$FDRthresh)
  out[["deDist"]] <- calcDist(deVS=out[["deVS"]])
  out[["deNeighb"]] <- calcNeighb(deVS=out[["deVS"]],
                                  deDist=out[["deDist"]])
  return(out)
}


#' Cluster-wise gene statistics
#'
#' Generates clusterwise gene statistics for use in
#' \code{\link{clusterWiseDEtest}}. You probably don't need to use this unless
#' you're trying to customize \code{\link{clusterWiseDEtest}}.
#'
#' @param nge The log-normalized gene expression matrix.
#'
#' @param cl The factor with cluster assignments per cell (column of nge).
#'
#' @param exponent The log base of your normalized input data. Seurat
#'   normalization uses the natural log (set this to exp(1)), while other
#'   normalization methods generally use log2 (set this to 2).
#'
#' @param pseudocount The pseudocount added to all log-normalized values in your
#'   input data. Most methods use a pseudocount of 1 to eliminate log(0) errors.
#'
#' @return The function returns a list of dataframes. Each list element contains
#'   a named list of clusters at that resolution. Each of those list elements
#'   contains a dataframe of three variables, where each sample is a gene.
#'   \code{DR} is the proportion of cells in the cluster in which that gene was
#'   detected. \code{MDTC} is mean normalized gene expression for that gene in
#'   only the cells in which it was detected (see \link{meanLogX} for mean
#'   calculation). \code{MTC} is the mean normalized gene expression for that
#'   gene in all cells of the cluster (see \link{meanLogX} for mean
#'   calculation).
#'
#' @seealso \code{\link{clusterWiseDEtest}} for the wrapper function to
#'   calculate all gene stats and differential expression tests for all input
#'   clusters.
#'
#' @export

calcCGS <- function(nge,cl,exponent,pseudocount) {
  print("-- Gene detection rate per cluster --")
  DR <- pbapply::pbapply(nge,1,function(X) tapply(X,cl,function(Y) sum(Y>0)/length(Y)))
  print("-- Mean detected gene expression per cluster --")
  MDTC <- pbapply::pbapply(nge,1,function(X) tapply(X,cl,function(Y) {
    temp <- meanLogX(Y[Y>0],ncell=ncol(nge),ex=exponent,pc=pseudocount)
    if (is.na(temp)) { temp <- 0 }
    return(temp)
  }))
  print("-- Mean gene expression per cluster --")
  MTC <- pbapply::pbapply(nge,1,function(X) 
    tapply(X,cl,function(Y) 
      meanLogX(Y,ncell=ncol(nge),ex=exponent,pc=pseudocount)))
  return(sapply(levels(cl),function(X) 
    data.frame(DR=DR[X,],MDTC=MDTC[X,],MTC=MTC[X,]),simplify=F))
}


#'Calculate logGER for deTissue calculation
#'
#'Calculates the log-ratios of gene expression for all genes in each one-vs-all
#'comparison of a cluster vs the rest of the data. This is used to determine the
#'genes used in deTissue calculations. You probably don't need to use this
#'unless you're trying to customize \code{\link{clusterWiseDEtest}}.
#'
#'@param nge The log-normalized gene expression matrix.
#'
#'@param cl The factor with cluster assignments per cell (column of nge).
#'
#'@param CGS The output from \code{\link{calcCGS}}.
#'
#'@param exponent The log base of your normalized input data. Seurat
#'  normalization uses the natural log (set this to exp(1)), while other
#'  normalization methods generally use log2 (set this to 2).
#'
#'@param pseudocount The pseudocount added to all log-normalized values in your
#'  input data. Most methods use a pseudocount of 1 to eliminate log(0) errors.
#'
#'@param logGERthresh Magnitude of gene expression ratio for a gene between
#'  clusters to use as filter for determining which genes to test for
#'  differential expression between clusters.
#'
#'@return The function returns a list where each list element is the log-ratios
#'  of gene expression when comparing each gene in a cluster to the rest of the
#'  cells as a whole in a one vs all comparison. These logGER tables are
#'  filtered to only include those gene that pass logGER threshold, and thus the
#'  names for each list entry correspond to the genes to test in
#'  \code{\link{calcDEvsRest}}.
#'
#'@seealso \code{\link{calcDEvsRest}} for the followup DE calculation, and
#'  \code{\link{clusterWiseDEtest}} for the wrapper function for all DE testing.
#'
#'@export

calcESvsRest <- function(nge,cl,CGS,exponent,pseudocount,logGERthresh) {
  deTes <- pbapply::pbsapply(levels(cl),function(i) 
    CGS[[i]]$MTC - apply(nge[,(cl != i | is.na(cl))],1,function(Y) 
      meanLogX(Y,ncell=ncol(nge),ex=exponent,pc=pseudocount)))
  return(sapply(colnames(deTes),function(X) 
    deTes[deTes[,X] > logGERthresh,X]))
}


#'Do deTissue calculation
#'
#'Calculates Wilcoxon rank-sum tests for all genes in each one-vs-all comparison
#'of a cluster vs the rest of the data. You probably don't need to use this
#'unless you're trying to customize \code{\link{clusterWiseDEtest}}.
#'
#'@param nge The log-normalized gene expression matrix.
#'
#'@param cl The factor with cluster assignments per cell (column of nge).
#'
#'@param deTes The output from \code{\link{calcESvsRest}}.
#'
#'@param exponent The log base of your normalized input data. Seurat
#'  normalization uses the natural log (set this to exp(1)), while other
#'  normalization methods generally use log2 (set this to 2).
#'
#'@param pseudocount The pseudocount added to all log-normalized values in your
#'  input data. Most methods use a pseudocount of 1 to eliminate log(0) errors.
#'
#'@param FDRthresh The false discovery rate to use as a threshold for
#'  determining statistical significance of differential expression calculated
#'  by the Wilcoxon rank-sum test.
#'
#'@return Differential testing results from Wilcoxon rank sum tests comparing a
#'  gene in each cluster to the rest of the cells as a whole in a one vs all
#'  comparison. The results are stored as a named list of dataframes. There is a
#'  list element for each cluster containing a dataframe of three variables,
#'  where each sample is a gene. \code{logGER} is the log gene expression ratio
#'  calculated by subtracting the mean expression of the gene (see
#'  \link{meanLogX} for mean calculation) in all other cells from the mean
#'  expression of the gene in this cluster. \code{pVal} is the p-value of the
#'  Wilcoxon rank sum test. \code{qVal} is the false discovery rate-corrected
#'  p-value of the test.
#'
#'@seealso \code{\link{calcESvsRest}} for the preceding calculation, and
#'  \code{\link{clusterWiseDEtest}} for the wrapper function for all DE testing.
#'
#'@export

calcDEvsRest <- function(nge,cl,deTes,exponent,pseudocount,FDRthresh) {
  deT_pVal <- pbapply::pbsapply(levels(cl),function(i)
    apply(nge[names(deTes[[i]]),],1,function(X) 
      suppressWarnings(wilcox.test(X[cl == i],X[cl != i])$p.value)),simplify=F)
  tempOut <- sapply(levels(cl),function(i) 
    data.frame(logGER=deTes[[i]],pVal=deT_pVal[[i]],
               qVal=p.adjust(deT_pVal[[i]],"fdr"))[order(deT_pVal[[i]]),],
    simplify=F)
  return(sapply(tempOut,function(X) X[X$qVal <= FDRthresh,],simplify=F))
}


#'Calculate genes used for deVS calculation
#'
#'Calculates the log-ratios of gene expression or difference in detection rate
#'for all genes in each of the potential combinations of clusters to compare.
#'This is used to determine the genes used in deVS calculations. You probably
#'don't need to use this unless you're trying to customize
#'\code{\link{clusterWiseDEtest}}.
#'
#'@param cl The factor with cluster assignments per cell (column of nge).
#'
#'@param CGS The output from \code{\link{calcCGS}}.
#'
#'@param threshType Filtering genes for use in differential expression testing
#'  can be done multiple ways. We use an expression ratio filter for comparing
#'  each cluster to the rest of the tissue as a whole, but find that difference
#'  in detection rates works better when comparing clusters to each other. You
#'  can set threshType to \code{"logGER"} to use a gene expression ratio for all
#'  gene filtering, or leave it as default (\code{"dDR"}) to use difference in
#'  detection rate as the thresholding method when comparing clusters to each
#'  other.
#'
#'@param dDRthresh Magnitude of detection rate difference of a gene between
#'  clusters to use as filter for determining which genes to test for
#'  differential expression between clusters.
#'
#'@param logGERthresh Magnitude of gene expression ratio for a gene between
#'  clusters to use as filter for determining which genes to test for
#'  differential expression between clusters.
#'
#'@return The function returns a list where each list element is a dataframe
#'  with effect size statistics (log-ratios of gene expression and difference in
#'  detection rate) when comparing each gene in a cluster to the rest of the
#'  cells as a whole in a one vs all comparison. These dataframes are filtered
#'  to only include those gene that pass the relevant threshold, and thus the
#'  rownames for each list entry correspond to the genes to test in
#'  \code{\link{calcDEvsCombn}}.
#'
#'@seealso \code{\link{calcDEvsCombn}} for the followup DE calculation, and
#'  \code{\link{clusterWiseDEtest}} for the wrapper function for all DE testing.
#'
#'@export

calcESvsCombn <- function(cl,CGS,threshType,logGERthresh,dDRthresh) {
  combos <- combn(levels(cl),2)
  colnames(combos) <- apply(combos,2,function(X) paste(X,collapse="-"))
  deMes <- apply(combos,2,function(i) 
    data.frame(logGER=CGS[[i[1]]]$MTC - CGS[[i[2]]]$MTC,
               dDR=CGS[[i[1]]]$DR - CGS[[i[2]]]$DR))
  for (i in names(deMes)) {
    rownames(deMes[[i]]) <- rownames(CGS[[1]])
  }
  return(sapply(names(deMes),function(X) 
    deMes[[X]][abs(deMes[[X]][,threshType]) > switch(threshType,
                                                     dDR=dDRthresh,
                                                     logGER=logGERthresh),],
    simplify=F))
}


#'Do deVS calculation
#'
#'Calculates Wilcoxon rank-sum tests for all genes in each of the potential
#'combinations of clusters to compare. You probably don't need to use this
#'unless you're trying to customize \code{\link{clusterWiseDEtest}}.
#'
#'@param nge The log-normalized gene expression matrix.
#'
#'@param cl The factor with cluster assignments per cell (column of nge).
#'
#'@param deMes The output from \code{\link{deMarkerGenesUsed}}.
#'
#'@param threshType Filtering genes for use in differential expression testing
#'  can be done multiple ways. We use an expression ratio filter for comparing
#'  each cluster to the rest of the tissue as a whole, but find that difference
#'  in detection rates works better when comparing clusters to each other. You
#'  can set threshType to \code{"logGER"} to use a gene expression ratio for all
#'  gene filtering, or leave it as default (\code{"dDR"}) to use difference in
#'  detection rate as the thresholding method when comparing clusters to each
#'  other.
#'
#'@param FDRthresh The false discovery rate to use as a threshold for
#'  determining statistical significance of differential expression calculated
#'  by the Wilcoxon rank-sum test.
#'
#'@return Differential testing results from Wilcoxon rank sum tests comparing a
#'  gene in each cluster to that gene in every other cluster in a series of
#'  tests. The results are stored as a nested list of dataframes. Each list
#'  element contains a named list of clusters (cluster A). Each of those lists
#'  contains a named list of all the other clusters (cluster B). Each of those
#'  list elements contains a dataframe of four variables, where each sample is a
#'  gene. \code{dDR} is the difference in detection rate of that gene between
#'  the two clusters (DR[A] - DR[B]). \code{logGER} is the log gene expression
#'  ratio calculated by taking the difference in mean expression of the gene
#'  (see \link{meanLogX} for mean calculation) between the two clusters (MTC[A]
#'  - MTC[B]). \code{pVal} is the p-value of the Wilcoxon rank sum test.
#'  \code{qVal} is the false discovery rate-corrected p-value of the test.
#'
#'
#'@seealso \code{\link{deMarkerGenesUsed}} for the preceding calculation, and
#'  \code{\link{clusterWiseDEtest}} for the wrapper function for all DE testing.
#'
#'@export

calcDEvsCombn <- function(nge,cl,deMes,threshType,FDRthresh) {
  combosL <- strsplit(names(deMes),"-")
  deM_pVal <- pbapply::pbsapply(seq_along(combosL),function(i)
    apply(nge[rownames(deMes[[i]]),],1,function(X) 
      suppressWarnings(wilcox.test(X[cl == combosL[[i]][1]],
                                   X[cl == combosL[[i]][2]])$p.value)),simplify=F)
  for (i in seq_along(deMes)) {
    deMes[[i]]$pVal <- deM_pVal[[i]]
    deMes[[i]]$qVal <- p.adjust(deM_pVal[[i]],"fdr")
    deMes[[i]] <- deMes[[i]][order(deMes[[i]]$qVal),]
  } 
  return(sapply(levels(cl),function(i) {
    temp <- list()
    for (X in seq_along(combosL)) {
      if (! i %in% combosL[[X]]) {
        next
      } else if (which(combosL[[X]] == i) == 1) {
        temp[[combosL[[X]][2]]] <- deMes[[X]][deMes[[X]][,threshType] > 0 & 
                                                   deMes[[X]]$qVal <= FDRthresh,]
      } else if (which(combosL[[X]] == i) == 2) {
        temp[[combosL[[X]][1]]] <- deMes[[X]][deMes[[X]][,threshType] < 0 &
                                                   deMes[[X]]$qVal <= FDRthresh,]
        temp[[combosL[[X]][1]]]$dDR <- temp[[combosL[[X]][1]]]$dDR * -1
        temp[[combosL[[X]][1]]]$logGER <- temp[[combosL[[X]][1]]]$logGER * -1
      }
    }
    return(temp)
  },simplify=F))
}


#'Do deMarker calculation
#'
#'Identifies marker genes for each cluster from \code{\link{calcDEvsCombn}} output.
#'You probably don't need to use this unless you're trying to customize
#'\code{\link{clusterWiseDEtest}}.
#'
#'@param deVS The output from \code{\link{calcDEvsCombn}}.
#'
#'@param FDRthresh The false discovery rate to use as a threshold for
#'  determining statistical significance of differential expression calculated
#'  by the Wilcoxon rank-sum test.
#'
#'@return Differential testing results from Wilcoxon rank sum tests comparing a
#'  gene in each cluster to that gene in every other cluster in a series of
#'  tests, and filtering for only those genes that show significant positive
#'  differential expression versus all other clusters. The results are stored as
#'  a list of dataframes. Each list element contains a named list of clusters
#'  (cluster A). Each of those list elements contains a dataframe where
#'  variables represent comparisons to all the other clusters and each sample is
#'  a gene. For each other cluster (cluster B), there are three variables, named
#'  as follows: \code{vs.B.dDR} is the difference in detection rate of that gene
#'  between the two clusters (DR[A] - DR[B]). \code{vs.B.logGER} is the log gene
#'  expression ratio calculated by taking the difference in mean expression of
#'  the gene (see \link{meanLogX} for mean calculation) between the two clusters
#'  (MTC[A] - MTC[B]). \code{vs.B.qVal} is the false discovery rate-corrected
#'  p-value of the Wilcoxon rank sum test.
#'
#'@seealso \code{\link{calcDEvsCombn}} for the preceding calculation, and
#'  \code{\link{clusterWiseDEtest}} for the wrapper function for all DE testing.
#'
#'@export

calcMarker <- function(deVS,FDRthresh) {
  return(sapply(deVS,function(X) {
    tempQval <- tapply(p.adjust(do.call(rbind,X)$pVal,"fdr"),
                       rep(names(X),sapply(X,nrow)),c,simplify=F)
    for (l in names(X)) {
      if (l %in% names(tempQval)) {
        X[[l]]$qVal <- tempQval[[l]]
        X[[l]] <- X[[l]][X[[l]]$qVal <= FDRthresh,]
      }
    }
    markerGenes <- Reduce(intersect,lapply(X,rownames))
    temp <- sapply(X,function(Y) Y[markerGenes,c("dDR","logGER","qVal")],simplify=F)
    names(temp) <- paste("vs",names(temp),sep=".")
    return(do.call(cbind,temp))
  },simplify=F))
}


#'Do deDist calculation
#'
#'Counts number of DE genes between each pair of clusters from
#'\code{\link{calcDEvsCombn}} output. You probably don't need to use this unless
#'you're trying to customize \code{\link{clusterWiseDEtest}}.
#'
#'@param deVS The output from \code{\link{calcDEvsCombn}}.
#'
#'@return A matrix of distances between clusters for each cluster resolution.
#'  Distances are calculated as number of differentially expressed genes between
#'  clusters. Interpretable by \code{\link{dist}} to generate a \code{dist}
#'  object.
#'
#'@seealso \code{\link{calcDEvsCombn}} for the preceding calculation, and
#'  \code{\link{clusterWiseDEtest}} for the wrapper function for all DE testing.
#'
#'@export

calcDist <- function(deVS) {
  sapply(names(deVS),function(X) sapply(names(deVS),function(Y) 
    if (X == Y) { return(NA) } else { min(nrow(deVS[[X]][[Y]]),nrow(deVS[[Y]][[X]])) }))
}


#'Do deNieghb calculation
#'
#'Identifies DE genes between nearest neighbour clusters using
#'\code{\link{calcDEvsCombn}} output for DE results and \code{\link{deDist}} output
#'to determine nearest neighbours. You probably don't need to use this unless
#'you're trying to customize \code{\link{clusterWiseDEtest}}.
#'
#'@param deVS The output from \code{\link{calcDEvsCombn}}.
#'
#'@param deDist The output from \code{\link{deDist}}.
#'
#'@return Differential testing results from Wilcoxon rank sum tests comparing a
#'  gene in each cluster to that gene in its nearest neighbouring cluster
#'  (calculated by number of differentially expressed genes), and filtering for
#'  only those genes that show significant positive differential expression
#'  versus all other clusters. The results are stored as a list of dataframes.
#'  Each list element contains a named list of clusters (cluster A). Each of
#'  those list elements contains a dataframe where variables represent the
#'  comparison to its nearest neighbouring cluster (cluster B) and each sample
#'  is a gene. There are three variables, named as follows: \code{vs.B.dDR} is
#'  the difference in detection rate of that gene between the two clusters
#'  (DR[A] - DR[B]). \code{vs.B.logGER} is the log gene expression ratio
#'  calculated by taking the difference in mean expression of the gene (see
#'  \link{meanLogX} for mean calculation) between the two clusters (MTC[A] -
#'  MTC[B]). \code{vs.B.qVal} is the false discovery rate-corrected p-value of
#'  the Wilcoxon rank sum test.
#'
#'@seealso \code{\link{calcDEvsCombn}} for the preceding calculation, and
#'  \code{\link{clusterWiseDEtest}} for the wrapper function for all DE testing.
#'
#'@export

calcNeighb <- function(deVS,deDist) {
  nb <- colnames(deDist)[apply(deDist,1,which.min)]
  names(nb) <- colnames(deDist)
  tempOut <- mapply(function(NB,VS) 
    VS[[NB]][,c("dDR","logGER","qVal")],
    NB=nb,VS=deVS,SIMPLIFY=F)
  for (i in names(tempOut)) {
    colnames(tempOut[[i]]) <- paste("vs",nb[i],colnames(tempOut[[i]]),sep=".")
  }
  return(tempOut)
}

