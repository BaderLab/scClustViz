#' @include sCVdataClass.R dataAccess.R
NULL


#' Prepare all cluster solutions for visualization with scClustViz
#'
#' An all-in-one function to prepare your data for viewing in the interactive
#' Shiny app. See example for the basic usage of scClustViz.
#'
#' This is a wrapper function for running \code{\link{CalcSCV}} over each
#' cluster resolution in the input, and outputs a list of \code{\link{sCVdata}}
#' objects that should be saved along with the input data. The resulting file is
#' ready to be read by \code{\link{runShiny}} for viewing. For each cluster
#' solution provided, this function calculates summary statistics per gene per
#' cluster, differential gene expression, and cluster separation metrics. This
#' may take a while to run, depending on the number of cluster solutions tested.
#' Use the \code{testAll} argument to prevent testing of overfitted cluster
#' solutions. To help track its progress, this function uses progress bars from
#' \code{pbapply}. To disable these, set
#' \code{\link[pbapply]{pboptions}(type="none")}. To re-enable, set
#' \code{\link[pbapply]{pboptions}(type="timer")}.
#'
#' @param inD The input dataset. An object of class \code{\link[Seurat]{seurat}}
#'   or \code{\link[SingleCellExperiment]{SingleCellExperiment}}. Other data
#'   classes are not currently supported.
#'   \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#'   for other data objects here!}
#' @param clusterDF A data frame of cluster assignments for all cells in the
#'   dataset. Variables (columns) are cluster solutions with different
#'   parameters, and rows should correspond to cells of the input gene
#'   expression matrix.
#' @param assayType Default = "" (for Seurat v1/2). A length-one character
#'   vector representing the assay slot in which the expression data is stored
#'   in the input object. This is not required for Seurat v1 or v2 objects. See
#'   \code{\link{getExpr}} for details.
#' @param DRforClust Default = "pca".A length-one character vector representing
#'   the dimensionality reduction method used as the input for clustering. This
#'   is commonly PCA, and should correspond to the slot name of the cell
#'   embedding in your input data - either the \code{type} argument in
#'   \code{\link[SingleCellExperiment]{reducedDim}(x,type)} or the
#'   \code{reduction.type} argument in
#'   \code{\link[Seurat]{GetDimReduction}(object,reduction.type)} (v2) or
#'   \code{reduction} in \code{\link[Seurat]{Embeddings}(object,reduction)}.
#' @param exponent Default = 2. A length-one numeric vector representing the
#'   base of the log-normalized gene expression data to be processed. Generally
#'   gene expression data is transformed into log2 space when normalizing (set
#'   this to 2), though \code{Seurat} uses the natural log (set this to exp(1)).
#' @param pseudocount Default = 1. A length-one numeric vector representing the
#'   pseudocount added to all log-normalized values in your input data. Most
#'   methods use a pseudocount of 1 to eliminate log(0) errors.
#' @param DRthresh Default = 0.1. A length-one numeric vector between 0 and 1
#'   representing the detection rate threshold for inclusion of a gene in the
#'   differential expression testing. A gene will be included if it is detected
#'   in at least this proportion of cells in at least one of the clusters being
#'   compared.
#' @param testAll Default = TRUE. Logical value indicating whether to test all
#'   cluster solutions (\code{TRUE}) or stop testing once a cluster solution has
#'   been found where there is no differentially expressed genes found between
#'   at least one pair of nearest neighbouring clusters (\code{FALSE}). If set
#'   to FALSE, this function will test cluster solutions in ascending order of
#'   number of clusters found. \emph{If set to (\code{FALSE}), only tested
#'   cluster solutions will appear in the scClustViz shiny app.}
#' @param FDRthresh Default = 0.05. A length-one numeric vector representing the
#'   targeted false discovery rate used to determine the number of
#'   differentially expressed genes between nearest neighbouring clusters,
#'   assuming \code{testAll} is set FALSE If \code{testAll} is TRUE, this
#'   argument is unused.
#' @param storeAllDE Default = TRUE. A logical vector of length 1 indicating
#'   whether to calculate and store effect size information for all genes in the
#'   comparison (TRUE), or just those passing the detection rate threshold for
#'   the Wilcoxon rank-sum test (FALSE). Setting this to FALSE will reduce the
#'   size of the output sCVdata object.
#' @param calcSil Default = TRUE. A logical vector of length 1. If TRUE,
#'   silhouette widths (a cluster cohesion/separation metric) will be calculated
#'   for all cells. This calculation is performed using the function
#'   \code{\link{CalcSilhouette}}, which is a wrapper to
#'   \code{\link[cluster]{silhouette}} with distance calculated using the same
#'   reduced dimensional cell embedding as was used for clustering, as indicated
#'   in the \code{DRforClust} argument. If the package \code{cluster} is not
#'   installed, this calculation is skipped.
#' @param calcDEvsRest Default = TRUE. A logical vector of length 1. If TRUE,
#'   differential expression tests will be performed comparing each cluster to
#'   the remaining cells in the data using a Wilcoxon rank-sum test and
#'   reporting false discovery rates. This calculation is performed using the
#'   function \code{\link{CalcDEvsRest}}. If set to FALSE, it is suggested that
#'   you perform DE testing on the same set of comparisons using a statistical
#'   method of your choice. This can be passed into your \code{sCVdata} objects
#'   in the list returned by \code{CalcAllSCV} using the function
#'   \code{\link{CalcDEvsRest}}. See function documentation for details.
#' @param calcDEcombn Default = TRUE.  A logical vector of length 1. If TRUE,
#'   differential expression tests will be performed comparing all pairwise
#'   combinations of clusters using a Wilcoxon rank-sum test and reporting false
#'   discovery rates. This calculation is performed using the function
#'   \code{\link{calcDEcombn}}. If set to FALSE, it is suggested that you
#'   perform DE testing on the same set of comparisons using a statistical
#'   method of your choice. This can be passed into your \code{sCVdata} objects
#'   in the list returned by \code{CalcAllSCV} using the function
#'   \code{\link{calcDEcombn}}. See function documentation for details.
#' @param UseBiocParallel Default = FALSE. Very experimental implementation of
#'   BiocParallel for calculations. Not recommended.
#'
#' @return The function returns a list containing \code{\link{sCVdata}} objects
#'   for each cluster resolution (sample) in the \code{clusterDF} data frame.
#'   The output object and the \code{inD} object should be saved as an
#'   \code{.RData} file. That file is the input for \code{\link{runShiny}}, the
#'   scClustViz Shiny interaction visualization app. See example. For details of
#'   calculations performed / stored by this function, see
#'   \code{\link{sCVdata}}.
#'
#' @examples
#' \dontrun{
#' your_cluster_columns <- grepl("res[.0-9]+$",
#'                               names(getMD(your_scRNAseq_data_object)))
#' # ^ Finds the cluster columns of the metadata in a Seurat object.
#'
#' your_cluster_results <- getMD(your_scRNAseq_data_object)[your_cluster_columns]
#'
#' sCVdata_list <- CalcAllSCV(inD=your_scRNAseq_data_object,
#'                            clusterDF=your_cluster_results,
#'                            DRforClust="pca",
#'                            exponent=exp(1),
#'                            pseudocount=1,
#'                            DRthresh=0.1,
#'                            testAll=F,
#'                            FDRthresh=0.05,
#'                            calcSil=T,
#'                            calcDEvsRest=T,
#'                            calcDEcombn=T)
#'
#' save(your_scRNAseq_data_object,sCVdata_list,
#'      file="for_scClustViz.RData")
#'
#' runShiny(filePath="for_scClustViz.RData")
#' # ^ see ?runShiny for detailed argument list
#' }
#'
#' @seealso \code{\link{sCVdata}} for information on the output data class.
#'   \code{\link{CalcSCV}} to generate an \code{sCVdata} object for a single
#'   cluster solution. \code{\link{runShiny}} starts the interactive Shiny GUI
#'   to view the results of this testing.
#'
#' @export

CalcAllSCV <- function(inD,
                       clusterDF,
                       assayType="",
                       DRforClust="pca",
                       exponent=2,
                       pseudocount=1,
                       DRthresh=0.1,
                       testAll=TRUE,
                       FDRthresh=0.05,
                       storeAllDE=T,
                       calcSil=T,
                       calcDEvsRest=T,
                       calcDEcombn=T,
                       UseBiocParallel=F) {
  if (!is(inD)[1] %in% findMethodSignatures(getExpr)) {
    stop(paste(
      paste0("Input data object must be one of: ",
             paste(findMethodSignatures(getExpr),collapse=", "),
             "."),
      paste("Other input objects are not supported at this time,",
            "but please let me know what object class"),
      paste("you'd like supported at",
            "https://github.com/BaderLab/scClustViz/issues, thanks!"),
      sep="\n  "))
  }
  if (is.null(colnames(getExpr(inD,assayType))) | is.null(rownames(getExpr(inD,assayType)))) {
    stop("Gene expression matrix returned by 'getExpr(inD,assayType)' is missing col/rownames.")
  }
  if (is.data.frame(clusterDF)) {
    if (nrow(clusterDF) != ncol(getExpr(inD,assayType))) {
      stop(paste("clusterDF must be a data frame where each variable (column) is the cluster",
                 "assignments for all cells (rows) from each cluster solution tested.",sep="\n  "))
    }
  } else if (length(clusterDF) == ncol(getExpr(inD,assayType))) {
    clusterDF <- as.data.frame(clusterDF)
    colnames(clusterDF) <- "Clust"
  } else {
    stop(paste("clusterDF must be a data frame where each variable (column) is the cluster",
               "assignments for all cells (rows) from each cluster solution tested.",sep="\n  "))
  }
  if (!identical(rownames(clusterDF),colnames(getExpr(inD,assayType)))) {
    rownames(clusterDF) <- colnames(getExpr(inD,assayType))
  }
  
  # If testAll == F, cluster solutions are sorted in ascending order of number
  # of clusters found.
  if (!testAll) {
    message(paste("  Testing cluster solutions in ascending order of number of clusters found.",
                  "Testing will stop after finding a solution with 0 differentially expressed",
                  "genes between nearest neighbouring clusters, and the resulting list of",
                  "sCVdata objects will be in ascending order of number of clusters found.",
                  sep="\n  "))
    sortedClusts <- order(sapply(clusterDF,function(X) length(unique(X))))
    clusterDF <- clusterDF[sortedClusts]
  }
  
  # This loop iterates through every cluster solution, and does DE testing
  # between clusters to generate the DE metrics for assessing your clusters.
  # This takes some time. If testAll == FALSE, it will stop once no
  # significantly differentially-expressed genes are detected between nearest
  # neighbouring clusters.
  outList <- list()
  for (X in names(clusterDF)) {
    message(paste(" ",
                  "--------------------------------------",
                  "--------------------------------------",
                  paste("Processing cluster solution:",X),
                  "--------------------------------------",sep="\n"))
    temp <- clusterDF[[X]]
    names(temp) <- rownames(clusterDF)
    outList[[X]] <- CalcSCV(inD=inD,
                            cl=temp,
                            assayType=assayType,
                            DRforClust=DRforClust,
                            exponent=exponent,
                            pseudocount=pseudocount,
                            DRthresh=DRthresh,
                            storeAllDE=storeAllDE,
                            calcSil=calcSil,
                            calcDEvsRest=calcDEvsRest,
                            calcDEcombn=calcDEcombn,
                            UseBiocParallel=UseBiocParallel)
    if (!testAll) {
      if (min(sapply(DEneighb(outList[[X]],FDRthresh),nrow)) < 1) { break }
    }
  }
  
  return(outList)
}


#' Create sCVdata object with calculation results for a cluster
#'
#' Creates and populates a new \code{\link{sCVdata}} object for a single
#' clustering result on the input data. This is the building block of the input
#' to \code{\link{runShiny}}, the scClustViz interactive visualization of
#' scRNAseq data.
#'
#' By default, \code{CalcSCV} populates all slots of the object, calculating
#' cluster separation metrics and differential gene expression results. At a
#' minimum, it calculates the basic gene statistics required by scClustViz for
#' visualization. This can take some time. To help track its progress, this
#' function uses progress bars from \code{pbapply}. To disable these, set
#' \code{\link[pbapply]{pboptions}(type="none")}. To re-enable, set
#' \code{\link[pbapply]{pboptions}(type="timer")}. To view the results using
#' \code{\link{runShiny}}, the resulting \code{\link{sCVdata}} object(s) must be
#' stored as a named list of cluster solutions and saved to an \code{.RData}
#' file along with the input data object. See example for details.
#'
#' @param inD The input dataset. An object of class \code{\link[Seurat]{seurat}}
#'   or \code{\link[SingleCellExperiment]{SingleCellExperiment}}. Other data
#'   classes are not currently supported.
#'   \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#'   for other data objects here!}
#' @param cl a factor where each value is the cluster assignment for a cell
#'   (column) in the input gene expression matrix.
#' @param assayType Default = "" (for Seurat v1/2). A length-one character
#'   vector representing the assay slot in which the expression data is stored
#'   in the input object. This is not required for Seurat v1 or v2 objects. See
#'   \code{\link{getExpr}} for details.
#' @param DRforClust Default = "pca".A length-one character vector representing
#'   the dimensionality reduction method used as the input for clustering. This
#'   is commonly PCA, and should correspond to the slot name of the cell
#'   embedding in your input data - either the \code{type} argument in
#'   \code{\link[SingleCellExperiment]{reducedDim}(x,type)} or the
#'   \code{reduction.type} argument in
#'   \code{\link[Seurat]{GetDimReduction}(object,reduction.type)} (v2) or
#'   \code{reduction} in \code{\link[Seurat]{Embeddings}(object,reduction)}.
#' @param exponent Default = 2. A length-one numeric vector representing the
#'   base of the log-normalized gene expression data to be processed. Generally
#'   gene expression data is transformed into log2 space when normalizing (set
#'   this to 2), though \code{Seurat} uses the natural log (set this to exp(1)).
#' @param pseudocount Default = 1. A length-one numeric vector representing the
#'   pseudocount added to all log-normalized values in your input data. Most
#'   methods use a pseudocount of 1 to eliminate log(0) errors.
#' @param DRthresh Default = 0.1. A length-one numeric vector between 0 and 1
#'   representing the detection rate threshold for inclusion of a gene in the
#'   differential expression testing. A gene will be included if it is detected
#'   in at least this proportion of cells in at least one of the clusters being
#'   compared.
#' @param storeAllDE Default = TRUE. A logical vector of length 1 indicating
#'   whether to calculate and store effect size information for all genes in the
#'   comparison (TRUE), or just those passing the detection rate threshold for
#'   the Wilcoxon rank-sum test (FALSE). Setting this to FALSE will reduce the
#'   size of the output sCVdata object.
#' @param calcSil Default = TRUE. A logical vector of length 1. If TRUE,
#'   silhouette widths (a cluster cohesion/separation metric) will be calculated
#'   for all cells. This calculation is performed using the function
#'   \code{\link{CalcSilhouette}}, which is a wrapper to
#'   \code{\link[cluster]{silhouette}} with distance calculated using the same
#'   reduced dimensional cell embedding as was used for clustering, as indicated
#'   in the \code{DRforClust} argument. If the package \code{cluster} is not
#'   installed, this calculation is skipped.
#' @param calcDEvsRest Default = TRUE. A logical vector of length 1. If TRUE,
#'   differential expression tests will be performed comparing each cluster to
#'   the remaining cells in the data using a Wilcoxon rank-sum test and
#'   reporting false discovery rates. This calculation is performed using the
#'   function \code{\link{CalcDEvsRest}}. If set to FALSE, it is suggested that
#'   you perform DE testing on the same set of comparisons using a statistical
#'   method of your choice. This can be passed into your \code{sCVdata} objects
#'   in the list returned by \code{CalcAllSCV} using the function
#'   \code{\link{CalcDEvsRest}}. See function documentation for details.
#' @param calcDEcombn Default = TRUE.  A logical vector of length 1. If TRUE,
#'   differential expression tests will be performed comparing all pairwise
#'   combinations of clusters using a Wilcoxon rank-sum test and reporting false
#'   discovery rates. This calculation is performed using the function
#'   \code{\link{calcDEcombn}}. If set to FALSE, it is suggested that you
#'   perform DE testing on the same set of comparisons using a statistical
#'   method of your choice. This can be passed into your \code{sCVdata} objects
#'   in the list returned by \code{CalcAllSCV} using the function
#'   \code{\link{calcDEcombn}}. See function documentation for details.
#' @param UseBiocParallel Default = FALSE. Very experimental implementation of
#'   BiocParallel for calculations. Not recommended.
#'
#' @return The function returns an \code{\link{sCVdata}} object with all slots
#'   populated by default, and at least the \code{Clusters},
#'   \code{ClustGeneStats}, and \code{params} slots populated. The resulting
#'   object can be added to a list as part of the input to
#'   \code{\link{runShiny}} for visualization of the cluster results.
#'
#' @examples
#' \dontrun{
#' ## This example shows integration of scClustViz with Seurat clustering ##
#'
#' DE_bw_clust <- TRUE
#' seurat_resolution <- 0
#' sCVdata_list <- list()
#'
#' while(DE_bw_clust) {
#'   seurat_resolution <- seurat_resolution + 0.2
#'   # ^ iteratively incrementing resolution parameter
#'
#'   your_seurat_obj <- Seurat::FindClusters(your_seurat_obj,
#'                                           resolution=seurat_resolution)
#'
#'   curr_sCVdata <- CalcSCV(inD=your_seurat_obj,
#'                           cl=your_seurat_obj@ident,
#'                           DRforClust="pca",
#'                           exponent=exp(1),
#'                           pseudocount=1,
#'                           DRthresh=0.1,
#'                           calcSil=T,
#'                           calcDEvsRest=T,
#'                           calcDEcombn=T)
#'
#'   DE_bw_NN <- sapply(DEneighb(curr_sCVdata,0.05),length)
#'   # ^ counts # of DE genes between neighbouring clusters at 5% FDR
#'
#'   if (min(DE_bw_NN) < 1) { DE_bw_clust <- FALSE }
#'   # ^ If no DE genes between nearest neighbours, don't loop again.
#'
#'   sCVdata_list[[paste0("res.",seurat_resolution)]] <- curr_sCVdata
#' }
#'
#' save(your_seurat_obj,sCVdata_list,
#'      file="for_scClustViz.RData")
#'
#' runShiny(filePath="for_scClustViz.RData")
#' # ^ see ?runShiny for detailed argument list
#' }
#'
#' @seealso \code{\link{sCVdata}} for information on the output data class.
#'   \code{\link{CalcAllSCV}} to generate a list of \code{sCVdata} objects for
#'   all cluster solutions in the input data. \code{\link{runShiny}} starts the
#'   interactive Shiny GUI to view the results of this testing.
#'
#' @export

CalcSCV <- function(inD,
                    cl,
                    assayType="",
                    DRforClust="pca",
                    exponent=2,
                    pseudocount=1,
                    DRthresh=0.1,
                    storeAllDE=T,
                    calcSil=T,
                    calcDEvsRest=T,
                    calcDEcombn=T,
                    UseBiocParallel=F) {
  if (!is(inD)[1] %in% findMethodSignatures(getExpr)) {
    stop(paste(
      paste0("Input data object must be one of: ",
             paste(findMethodSignatures(getExpr),collapse=", "),
             "."),
      paste("Other input objects are not supported at this time,",
            "but please let me know what object class"),
      paste("you'd like supported at",
            "https://github.com/BaderLab/scClustViz/issues, thanks!"),
      sep="\n  "))
  }
  if (is.null(colnames(getExpr(inD,assayType))) | is.null(rownames(getExpr(inD,assayType)))) {
    stop("Gene expression matrix returned by 'getExpr(inD,assayType)' is missing col/rownames.")
  }
  if (length(cl) != ncol(getExpr(inD,assayType))) {
    stop(paste("cl must be a factor where each value is the cluster assignment",
               "for a cell (column) in the input gene expression matrix.",
               sep="\n  "))
  }
  if (is.character(cl)) {
    cl <- as.factor(cl)
  }
  if (!all(names(cl) == colnames(getExpr(inD,assayType))) | is.null(names(cl))) {
    names(cl) <- colnames(getExpr(inD,assayType))
  }
  if (any(grepl("-",levels(cl)))) {
    stop("Cluster names (levels in 'cl') cannot contain '-'.")
  }
  
  out <- sCVdata(Clusters=cl,
                 params=sCVparams(assayType=assayType,
                                  DRforClust=DRforClust,
                                  exponent=exponent,
                                  pseudocount=pseudocount,
                                  DRthresh=DRthresh))
  if (calcSil) { #Doing this first in case DRforClust is set incorrectly.
    if (require(cluster)) {
      Silhouette(out) <- CalcSilhouette(out,inD)
    } else {
      warning(paste("  Silhouette could not be calculated because package 'cluster'",
                    "is not installed. Try 'install.packages(cluster)', then run",
                    "'CalcSilhouette()' for the sCVdata object at this resolution.",sep="\n  "))
    }
  }
  
  #this is not optional, since everything depends on it.
  ClustGeneStats(out) <- CalcCGS(out,inD,UseBiocParallel) 
  
  if (calcDEvsRest) {
    DEvsRest(out) <- CalcDEvsRest(out,inD,storeAllDE,UseBiocParallel)
  }
  
  if (calcDEcombn) {
    DEcombn(out) <- CalcDEcombn(out,inD,storeAllDE,UseBiocParallel)
  }
  return(out)
}


# CalcCGS ----

#' Internal fx for cluster-wise gene statistics
#'
#' Internal function. See \code{\link{CalcCGS}}.
#'
#' @param nge The log-normalized gene expression matrix.
#' @param cl The factor with cluster assignments per cell (column of nge).
#' @param exponent The log base of your normalized input data. Seurat
#'   normalization uses the natural log (set this to exp(1)), while other
#'   normalization methods generally use log2 (set this to 2).
#' @param pseudocount The pseudocount added to all log-normalized values in your
#'   input data. Most methods use a pseudocount of 1 to eliminate log(0) errors.
#' 
#' @return The function returns a list of dataframes. Each list element contains
#'   a named list of clusters at that resolution. Each of those list elements
#'   contains a dataframe of three variables, where each sample is a gene.
#'   \code{DR} is the proportion of cells in the cluster in which that gene was
#'   detected. \code{MDGE} is mean normalized gene expression for that gene in
#'   only the cells in which it was detected (see \code{\link{meanLogX}} for
#'   mean calculation). \code{MGE} is the mean normalized gene expression for
#'   that gene in all cells of the cluster (see \code{\link{meanLogX}} for mean
#'   calculation).

fx_calcCGS <- function(nge,cl,exponent,pseudocount) {
  message("-- Calculating gene detection rate per cluster --")
  DR <- pbapply::pbsapply(sapply(levels(cl),function(i) nge[,cl %in% i,drop=F],simplify=F),
                          function(X) apply(X,1,function(Y) sum(Y > 0)/length(Y)),simplify=F)
  
  message("-- Calculating mean detected gene expression per cluster --")
  MDGE <- pbapply::pbsapply(sapply(levels(cl),function(i) nge[,cl %in% i,drop=F],simplify=F),
                            function(X) apply(X,1,function(Y) {
                              temp <- meanLogX(Y[Y > 0],
                                               ncell=ncol(nge),
                                               ex=exponent,
                                               pc=pseudocount)
                              if (is.na(temp)) { temp <- 0 }
                              return(temp)
                            }),simplify=F)
  
  message("-- Calculating mean gene expression per cluster --")
  MGE <- pbapply::pbsapply(sapply(levels(cl),function(i) nge[,cl %in% i,drop=F],simplify=F),
                           function(X) apply(X,1,function(Y)
                             meanLogX(Y,
                                      ncell=ncol(nge),
                                      ex=exponent,
                                      pc=pseudocount)),simplify=F)
  
  return(sapply(levels(cl),function(X) 
    data.frame(DR=DR[[X]],MDGE=MDGE[[X]],MGE=MGE[[X]]),simplify=F))
}

#' Internal fx for cluster-wise gene statistics using BiocParallel
#'
#' Internal function. See \code{\link{CalcCGS}}.
#'
#' @param nge The log-normalized gene expression matrix.
#' @param cl The factor with cluster assignments per cell (column of nge).
#' @param exponent The log base of your normalized input data. Seurat
#'   normalization uses the natural log (set this to exp(1)), while other
#'   normalization methods generally use log2 (set this to 2).
#' @param pseudocount The pseudocount added to all log-normalized values in your
#'   input data. Most methods use a pseudocount of 1 to eliminate log(0) errors.
#' 
#' @return The function returns a list of dataframes. Each list element contains
#'   a named list of clusters at that resolution. Each of those list elements
#'   contains a dataframe of three variables, where each sample is a gene.
#'   \code{DR} is the proportion of cells in the cluster in which that gene was
#'   detected. \code{MDGE} is mean normalized gene expression for that gene in
#'   only the cells in which it was detected (see \code{\link{meanLogX}} for
#'   mean calculation). \code{MGE} is the mean normalized gene expression for
#'   that gene in all cells of the cluster (see \code{\link{meanLogX}} for mean
#'   calculation).

fx_calcCGS_BP <- function(nge,cl,exponent,pseudocount) {
  message("-- Calculating gene detection rate per cluster --")
  DR <- BiocParallel::bplapply(sapply(levels(cl),function(i) nge[,cl %in% i],simplify=F),
                               function(X) apply(X,1,function(Y) sum(Y>0)/length(Y)))
  names(DR) <- levels(cl)

  message("-- Calculating mean detected gene expression per cluster --")
  MDGE <- BiocParallel::bplapply(sapply(levels(cl),function(i) nge[,cl %in% i],simplify=F),
                                 function(X) apply(X,1,function(Y) {
                                   temp <- meanLogX(Y[Y>0],
                                                    ncell=ncol(nge),
                                                    ex=exponent,
                                                    pc=pseudocount)
                                   if (is.na(temp)) { temp <- 0 }
                                   return(temp)
                                 }))
  names(MDGE) <- levels(cl)

  message("-- Calculating mean gene expression per cluster --")
  MGE <- BiocParallel::bplapply(sapply(levels(cl),function(i) nge[,cl %in% i],simplify=F),
                                function(X) apply(X,1,function(Y)
                                  meanLogX(Y,
                                           ncell=ncol(nge),
                                           ex=exponent,
                                           pc=pseudocount)))
  names(MGE) <- levels(cl)

  return(sapply(levels(cl),function(X) 
    data.frame(DR=DR[[X]],MDGE=MDGE[[X]],MGE=MGE[[X]]),simplify=F))
}


#' Calculate cluster-wise gene statistics for sCVdata
#'
#' Calculates gene summary statistics per cluster for the clusters in an sCVdata
#' object, using the gene expression matrix from the input data object. This is
#' called by \code{\link{CalcSCV}} and you shouldn't need to call it on its own.
#'
#' To help track its progress, this function uses progress bars from
#' \code{pbapply}. To disable these, set
#' \code{\link[pbapply]{pboptions}(type="none")}. To re-enable, set
#' \code{\link[pbapply]{pboptions}(type="timer")}.
#'
#' @param sCVd An sCVdata object.
#' @param inD The input dataset. An object of class \code{\link[Seurat]{seurat}}
#'   or \code{\link[SingleCellExperiment]{SingleCellExperiment}}. Other data
#'   classes are not currently supported.
#'   \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#'   for other data objects here!}
#' @param UseBiocParallel Default = FALSE. Very experimental implementation of
#'   BiocParallel for calculations. Not recommended.
#'
#' @return The function returns a list of dataframes. Each list element contains
#'   a named list of clusters at that resolution. Each of those list elements
#'   contains a dataframe where each sample is a gene, containing the following
#'   variables: \code{DR} is the proportion of cells in the cluster in which
#'   that gene was detected. \code{MDGE} is mean normalized gene expression for
#'   that gene in only the cells in which it was detected (see
#'   \code{\link{meanLogX}} for mean calculation). \code{MGE} is the mean
#'   normalized gene expression for that gene in all cells of the cluster (see
#'   \code{\link{meanLogX}} for mean calculation).
#'
#' @seealso \code{\link{CalcSCV}} for wrapper function to calculate all
#'   statistics for an sCVdata object, and \code{\link{fx_calcCGS}} for the
#'   internal function performing the calculations.
#'
#' @examples
#' \dontrun{
#' ClustGeneStats(your_sCV_obj) <- CalcCGS(sCVd=your_sCV_obj,
#'                                         inD=your_scRNAseq_data_object)
#' }
#' 
#' @name CalcCGS
#'
#' @export
#' 

setGeneric("CalcCGS",function(sCVd,inD,UseBiocParallel) standardGeneric("CalcCGS"))


#' @describeIn CalcCGS Calculate cluster-wise gene stats for sCVdata
#' @export

setMethod("CalcCGS",signature("sCVdata"),
          function(sCVd,inD,UseBiocParallel=FALSE) {
            if (UseBiocParallel) {
              fx_calcCGS_BP(nge=getExpr(inD,Param(sCVd,"assayType")),
                            cl=Clusters(sCVd),
                            exponent=Param(sCVd,"exponent"),
                            pseudocount=Param(sCVd,"pseudocount"))
            } else {
              fx_calcCGS(nge=getExpr(inD,Param(sCVd,"assayType")),
                         cl=Clusters(sCVd),
                         exponent=Param(sCVd,"exponent"),
                         pseudocount=Param(sCVd,"pseudocount"))
            }
          })


# CalcDEvsRest ----

#' Internal fx to calculate logGER for DEvsRest calculation
#'
#' Internal function. See \code{\link{CalcDEvsRest}}.
#'
#' Calculates the log-ratios of gene expression for all genes in each one-vs-all
#' comparison of a cluster vs the rest of the data. This is used to determine
#' the genes used in DEvsRest calculations.
#'
#' @param nge The log-normalized gene expression matrix.
#' @param cl The factor with cluster assignments per cell (column of nge).
#' @param CGS The output from \code{\link{CalcCGS}}.
#' @param exponent The log base of your normalized input data. Seurat
#'   normalization uses the natural log (set this to exp(1)), while other
#'   normalization methods generally use log2 (set this to 2).
#' @param pseudocount The pseudocount added to all log-normalized values in your
#'   input data. Most methods use a pseudocount of 1 to eliminate log(0) errors.
#' @param DRthresh The threshold for minimum detection rate of a gene in the
#'   cluster for the gene to be considered in the following Wilcoxon rank-sum
#'   test.
#'
#' @return The function returns a list where each list element is the log-ratios
#'   of gene expression when comparing each gene in a cluster to the rest of the
#'   cells as a whole in a one vs all comparison. These logGER tables are
#'   filtered to only include those gene that pass logGER threshold, and thus
#'   the names for each list entry correspond to the genes to test in
#'   \code{\link{fx_calcDEvsRest}}.
#'   

fx_calcESvsRest <- function(nge,cl,CGS,exponent,pseudocount,DRthresh) {
  message("-- Calculating differential expression cluster vs rest effect size --")
  return(pbapply::pbsapply(levels(cl),function(i) {
    temp <- data.frame(overThreshold = CGS[[i]]$DR >= DRthresh,
                       logGER=NA,
                       Wstat=NA,
                       pVal=NA,
                       FDR=NA)
    rownames(temp) <- rownames(CGS[[i]])
    temp[temp$overThreshold,"logGER"] <- CGS[[i]][temp$overThreshold,"MGE"] - 
      apply(nge[temp$overThreshold,(!cl %in% i | is.na(cl))],1,function(Y) 
        meanLogX(Y,ncell=ncol(nge),ex=exponent,pc=pseudocount))
    return(temp)
  },simplify=F))
}


#' Internal fx to calculate logGER for DEvsRest calculation using BiocParallel
#'
#' Internal function. See \code{\link{CalcDEvsRest}}.
#'
#' Calculates the log-ratios of gene expression for all genes in each one-vs-all
#' comparison of a cluster vs the rest of the data. This is used to determine
#' the genes used in DEvsRest calculations.
#'
#' @param nge The log-normalized gene expression matrix.
#' @param cl The factor with cluster assignments per cell (column of nge).
#' @param CGS The output from \code{\link{CalcCGS}}.
#' @param exponent The log base of your normalized input data. Seurat
#'   normalization uses the natural log (set this to exp(1)), while other
#'   normalization methods generally use log2 (set this to 2).
#' @param pseudocount The pseudocount added to all log-normalized values in your
#'   input data. Most methods use a pseudocount of 1 to eliminate log(0) errors.
#' @param DRthresh The threshold for minimum detection rate of a gene in the
#'   cluster for the gene to be considered in the following Wilcoxon rank-sum
#'   test.
#'
#' @return The function returns a list where each list element is the log-ratios
#'   of gene expression when comparing each gene in a cluster to the rest of the
#'   cells as a whole in a one vs all comparison. These logGER tables are
#'   filtered to only include those gene that pass logGER threshold, and thus
#'   the names for each list entry correspond to the genes to test in
#'   \code{\link{fx_calcDEvsRest}}.
#'   

fx_calcESvsRest_BP <- function(nge,cl,CGS,exponent,pseudocount,DRthresh) {
  message("-- Calculating differential expression cluster vs rest effect size --")
  temp <- BiocParallel::bplapply(levels(cl),function(i) {
    temp <- data.frame(overThreshold = CGS[[i]]$DR >= DRthresh,
                       logGER=NA,
                       Wstat=NA,
                       pVal=NA,
                       FDR=NA)
    rownames(temp) <- rownames(CGS[[i]])
    temp[temp$overThreshold,"logGER"] <- CGS[[i]][temp$overThreshold,"MGE"] -
      apply(nge[temp$overThreshold,(!cl %in% i | is.na(cl))],1,function(Y)
        meanLogX(Y,ncell=ncol(nge),ex=exponent,pc=pseudocount))
    return(temp)
  })
  names(temp) <- levels(cl)
  return(temp)
}


#' Internal fx to perform one vs all DE testing
#'
#' Internal function. See \code{\link{CalcDEvsRest}}.
#'
#' Calculates Wilcoxon rank-sum tests for all genes in each one-vs-all
#' comparison of a cluster vs the rest of the data. You probably don't need to
#' use this unless you're trying to customize \code{\link{clusterWiseDEtest}}.
#'
#' @param nge The log-normalized gene expression matrix.
#' @param cl The factor with cluster assignments per cell (column of nge).
#' @param deTes The output from \code{\link{fx_calcESvsRest}}.
#'
#' @return Differential testing results from Wilcoxon rank sum tests comparing a
#'   gene in each cluster to the rest of the cells as a whole in a one vs all
#'   comparison. The results are stored as a named list of dataframes. There is
#'   a list element for each cluster containing a dataframe of three variables,
#'   where each sample is a gene. \code{logGER} is the log gene expression ratio
#'   calculated by subtracting the mean expression of the gene (see
#'   \link{meanLogX} for mean calculation) in all other cells from the mean
#'   expression of the gene in this cluster. \code{Wstat} and \code{pVal} are
#'   the test statistic and the p-value of the Wilcoxon rank sum test.
#'   \code{FDR} is the false discovery rate-corrected p-value of the test.
#'   

fx_calcDEvsRest <- function(nge,cl,deTes) {
  message("-- Testing differential expression cluster vs rest --")
  deT_pVal <- presto::wilcoxauc(X=nge,y=cl)
  for (i in names(deTes)) {
    tempRows <- deT_pVal$feature %in% rownames(deTes[[i]])[deTes[[i]]$overThreshold] & deT_pVal$group == i
    deTes[[i]][deT_pVal[tempRows,"feature"],"Wstat"] <- deT_pVal[tempRows,"statistic"]
    deTes[[i]][deT_pVal[tempRows,"feature"],"pVal"] <- deT_pVal[tempRows,"pval"]
    deTes[[i]][deT_pVal[tempRows,"feature"],"FDR"] <- p.adjust(deT_pVal[tempRows,"pval"],"fdr")
  } 
  return(deTes)
}


#' Internal fx to perform one vs all DE testing using BiocParallel
#'
#' Internal function. See \code{\link{CalcDEvsRest}}.
#'
#' Calculates Wilcoxon rank-sum tests for all genes in each one-vs-all
#' comparison of a cluster vs the rest of the data. You probably don't need to
#' use this unless you're trying to customize \code{\link{clusterWiseDEtest}}.
#'
#' @param nge The log-normalized gene expression matrix.
#' @param cl The factor with cluster assignments per cell (column of nge).
#' @param deTes The output from \code{\link{fx_calcESvsRest}}.
#'
#' @return Differential testing results from Wilcoxon rank sum tests comparing a
#'   gene in each cluster to the rest of the cells as a whole in a one vs all
#'   comparison. The results are stored as a named list of dataframes. There is
#'   a list element for each cluster containing a dataframe of three variables,
#'   where each sample is a gene. \code{logGER} is the log gene expression ratio
#'   calculated by subtracting the mean expression of the gene (see
#'   \link{meanLogX} for mean calculation) in all other cells from the mean
#'   expression of the gene in this cluster. \code{Wstat} and \code{pVal} are
#'   the test statistic and the p-value of the Wilcoxon rank sum test.
#'   \code{FDR} is the false discovery rate-corrected p-value of the test.
#'   

fx_calcDEvsRest_BP <- function(nge,cl,deTes) {
  message("-- Testing differential expression cluster vs rest --")
  deT_pVal <- BiocParallel::bplapply(levels(cl),function(i)
    apply(nge[rownames(deTes[[i]])[deTes[[i]]$overThreshold],],1,function(X)
      # suppressWarnings(wilcox.test(X[cl %in% i],X[!cl %in% i],alternative="greater")$p.value)
      suppressWarnings(unlist(wilcox.test(X[cl %in% i],X[!cl %in% i])[c("statistic","p.value")]))
    ))
  names(deT_pVal) <- levels(cl)
  for (i in names(deTes)) {
    deTes[[i]][colnames(deT_pVal[[i]]),"Wstat"] <- deT_pVal[[i]]["statistic.W",]
    deTes[[i]][colnames(deT_pVal[[i]]),"pVal"] <- deT_pVal[[i]]["p.value",]
    deTes[[i]][colnames(deT_pVal[[i]]),"FDR"] <- p.adjust(deT_pVal[[i]]["p.value",],"fdr")
  } 
  return(deTes)
}


#' Calculates one vs. all DE tests for sCVdata
#'
#' Performs differential gene expression tests for each cluster in an sCVdata
#' object, comparing the cells in the cluster to the remaining cells in the data
#' using the gene expression matrix of input data object. Alternatively, this
#' function can be skipped, and existing DE test results can be assigned
#' directly to the sCVdata object.
#'
#' This function performs Wilcoxon rank sum tests comparing gene expression
#' between each cluster and all other cells in the input data. Gene expression
#' ratio in log space (\code{logGER}) is reported for all genes in the
#' comparison. Genes are tested if they are detected in the cluster at a higher
#' proportion than \code{Param(sCVd,"DRthresh")}, and both unadjusted p-values
#' and false discovery rates are reported for all genes tested. To help track
#' its progress, this function uses progress bars from \code{pbapply}. To
#' disable these, set \code{\link[pbapply]{pboptions}(type="none")}. To
#' re-enable, set \code{\link[pbapply]{pboptions}(type="timer")}.
#'
#' If using existing DE test results, assign results of one vs. all tests for
#' every cluster in sCVdata to the \code{\link{DEvsRest}} slot of the
#' \code{\link{sCVdata}} object. See example and slot documentation.
#'
#' @param sCVd An sCVdata object.
#' @param inD The input dataset. An object of class \code{\link[Seurat]{seurat}}
#'   or \code{\link[SingleCellExperiment]{SingleCellExperiment}}. Other data
#'   classes are not currently supported.
#'   \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#'   for other data objects here!}
#' @param storeAllDE A logical vector of length 1 indicating whether to
#'   calculate and store effect size information for all genes in the comparison
#'   (TRUE), or just those passing the detection rate threshold for the Wilcoxon
#'   rank-sum test (FALSE). Setting this to FALSE will reduce the size of the
#'   output sCVdata object.
#' @param UseBiocParallel Default = FALSE. Very experimental implementation of
#'   BiocParallel for calculations. Not recommended.
#'
#' @return A named list of data frames, one entry for each level in
#'   \code{Clusters(sCVd)} (with corresponding name). Each entry is data frame
#'   containing gene differential expression stats when comparing the cells of
#'   that cluster to all other cells in the input data. Rows represent genes,
#'   and variables include \code{logGER} (an effect size measure: gene
#'   expression ratio in log space, often referred to as logFC) and \code{FDR}
#'   (significance measure: false discovery rate). Also included are
#'   \code{Wstat} and \code{pVal}, the test statistic and the p-value of the
#'   Wilcoxon rank sum test.
#'
#' @seealso \code{\link{CalcSCV}} for wrapper function to calculate all
#'   statistics for an sCVdata object,  and \code{\link{fx_calcESvsRest}} and
#'   \code{\link{fx_calcDEvsRest}} for the internal functions performing the
#'   calculations. Wilcox test is now powered by \code{\link[presto]{wilcoxauc}}
#'   for super speed.
#'
#' @examples
#' \dontrun{
#' ## Example using CalcDEvsRest ##
#' DEvsRest(your_sCV_obj) <- CalcDEvsRest(sCVd=your_sCV_obj,
#'                                        inD=your_scRNAseq_data_object)
#'
#'
#' ## Example using MAST results from Seurat to replace CalcDEvsRest ##
#' MAST_oneVSall <- FindAllMarkers(your_seurat_obj,
#'                                 logfc.threshold=0,
#'                                 min.pct=0.1,
#'                                 test.use="MAST",
#'                                 latent.vars="nUMI")
#' # ^ FindAllMarkers and CalcDEvsRest do equivalent comparisons 
#' 
#' names(MAST_oneVSall)[names(MAST_oneVSall) == "avg_logFC"] <- "logGER"
#' # ^ Effect size variable must be named 'logGER'
#' names(MAST_oneVSall)[names(MAST_oneVSall) == "p_val_adj"] <- "FDR"
#' # ^ Significance variable must be named 'FDR'
#' 
#' MAST_oneVSall_list <- sapply(levels(MAST_oneVSall$cluster),
#'                              function(X) {
#'                                temp <- MAST_oneVSall[MAST_oneVSall$cluster == X,]
#'                                rownames(temp) <- temp$gene
#'                                # ^ Rownames must be gene names.
#'                                return(temp)
#'                              },simplify=F)
#' # ^ Dataframe converted to list of dataframes per cluster
#' 
#' DEvsRest(your_sCV_obj) <- MAST_oneVSall_list
#' # ^ Slot MAST results into sCVdata object
#'
#' }
#'
#' @name CalcDEvsRest
#'
#' @export
#' 

setGeneric("CalcDEvsRest",function(sCVd,inD,storeAllDE,UseBiocParallel) 
  standardGeneric("CalcDEvsRest"))


#' @describeIn CalcDEvsRest Calculate one vs. all DE tests for sCVdata
#' @export

setMethod("CalcDEvsRest","sCVdata",
          function(sCVd,inD,storeAllDE=TRUE,UseBiocParallel=FALSE) {
            if (!is(inD)[1] %in% findMethodSignatures(getExpr)) {
              stop(paste("The input data object must be one of:",
                         paste(findMethodSignatures(getExpr),collapse=", "),
                         sep="\n  "))
            }
            if (length(levels(Clusters(sCVd))) <= 1) {
              stop("scClustViz can't calculate differential expression when there's only one cluster.")
            }
            
            if (UseBiocParallel) {
              deTes <- fx_calcESvsRest_BP(nge=getExpr(inD,Param(sCVd,"assayType")),
                                          cl=Clusters(sCVd),
                                          CGS=ClustGeneStats(sCVd),
                                          exponent=Param(sCVd,"exponent"),
                                          pseudocount=Param(sCVd,"pseudocount"),
                                          DRthresh=Param(sCVd,"DRthresh"))
            } else {
              deTes <- fx_calcESvsRest(nge=getExpr(inD,Param(sCVd,"assayType")),
                                       cl=Clusters(sCVd),
                                       CGS=ClustGeneStats(sCVd),
                                       exponent=Param(sCVd,"exponent"),
                                       pseudocount=Param(sCVd,"pseudocount"),
                                       DRthresh=Param(sCVd,"DRthresh"))
            }
            if (!storeAllDE) { 
              deTes <- sapply(deTes,function(X) X[X$overThreshold,],simplify=F) 
            }
            if (UseBiocParallel) {
              deTes <- fx_calcDEvsRest_BP(nge=getExpr(inD,Param(sCVd,"assayType")),
                                          cl=Clusters(sCVd),
                                          deTes=deTes)
            } else {
              deTes <- fx_calcDEvsRest(nge=getExpr(inD,Param(sCVd,"assayType")),
                                       cl=Clusters(sCVd),
                                       deTes=deTes)
            }
            return(deTes)
          })


# CalcDEcombn ----

#' Internal fx to calculate genes used for deVS calculation
#'
#' Internal function. See \code{\link{CalcDEcombn}}.
#'
#' Calculates the log-ratios of gene expression and difference in detection rate
#' for all genes in each of the potential combinations of clusters to compare.
#' This is used to determine the genes used in deVS calculations.
#'
#' @param cl The factor with cluster assignments per cell (column of nge).
#' @param CGS The output from \code{\link{fx_calcCGS}}.
#' @param DRthresh The threshold for minimum detection rate of a gene to be
#'   considered in the following Wilcoxon rank-sum test. Gene must pass
#'   threshold in at least one of the pair of clusters.
#'
#' @return The function returns a list where each list element is a dataframe
#'   with effect size statistics (log-ratios of gene expression and difference
#'   in detection rate) when comparing each gene in a cluster to the rest of the
#'   cells as a whole in a one vs all comparison. These dataframes are filtered
#'   to only include those gene that pass the relevant threshold, and thus the
#'   rownames for each list entry correspond to the genes to test in
#'   \code{\link{fx_calcDEcombn}}.
#'

fx_calcEScombn <- function(cl,CGS,DRthresh) {
  combos <- combn(levels(cl),2)
  colnames(combos) <- apply(combos,2,function(X) paste(X,collapse="-"))
  return(apply(combos,2,function(i) {
    temp <- data.frame(overThreshold=(CGS[[i[1]]]$DR >= DRthresh | CGS[[i[2]]]$DR >= DRthresh),
                       logGER=CGS[[i[1]]]$MGE - CGS[[i[2]]]$MGE,
                       dDR=CGS[[i[1]]]$DR - CGS[[i[2]]]$DR)
    rownames(temp) <- rownames(CGS[[i[1]]])
    return(temp)
  }))
}


#' Internal fx to calculate DE between combinations of clusters
#'
#' Internal function. See \code{\link{CalcDEcombn}}.
#'
#' Calculates Wilcoxon rank-sum tests for all genes in each of the potential
#' combinations of clusters to compare.
#'
#' @param nge The log-normalized gene expression matrix.
#' @param cl The factor with cluster assignments per cell (column of nge).
#' @param deMes The output from \code{\link{fx_calcEScombn}}.
#'
#' @return Differential testing results from Wilcoxon rank sum tests comparing a
#'   gene in each cluster to that gene in every other cluster in a series of
#'   tests. The results are stored as a nested list of dataframes. Each list
#'   element contains a named list of clusters (cluster A). Each of those lists
#'   contains a named list of all the other clusters (cluster B). Each of those
#'   list elements contains a dataframe of four variables, where each sample is
#'   a gene. \code{dDR} is the difference in detection rate of that gene between
#'   the two clusters (DR[A] - DR[B]). \code{logGER} is the log gene expression
#'   ratio calculated by taking the difference in mean expression of the gene
#'   (see \code{\link{meanLogX}} for mean calculation) between the two clusters
#'   (MGE[A] - MGE[B]). \code{Wstat} and \code{pVal} are the test statistic and
#'   the p-value of the Wilcoxon rank sum test. \code{FDR} is the false
#'   discovery rate-corrected p-value of the test.
#'   

fx_calcDEcombn <- function(nge,cl,deMes) {
  combosL <- strsplit(names(deMes),"-")
  message("-- Testing differential expression between clusters --")
  deM_pVal <- pbapply::pbsapply(combosL,function(G) {
    temp <- presto::wilcoxauc(X=nge,y=cl,groups_use=G)
    temp <- temp[temp$group == G[1],]
    return(temp)
  },simplify=F)
  names(deM_pVal) <- names(deMes)
  for (i in names(deMes)) {
    tempRows <- deM_pVal[[i]]$feature %in% rownames(deMes[[i]])[deMes[[i]]$overThreshold]
    deMes[[i]][deM_pVal[[i]][tempRows,"feature"],"Wstat"] <- deM_pVal[[i]][tempRows,"statistic"]
    deMes[[i]][deM_pVal[[i]][tempRows,"feature"],"pVal"] <- deM_pVal[[i]][tempRows,"pval"]
    deMes[[i]][deM_pVal[[i]][tempRows,"feature"],"FDR"] <- p.adjust(deM_pVal[[i]][tempRows,"pval"],"fdr")
  } 
  return(deMes)
}


#' Internal fx to calculate DE between combinations of clusters using BiocParallel
#'
#' Internal function. See \code{\link{CalcDEcombn}}.
#'
#' Calculates Wilcoxon rank-sum tests for all genes in each of the potential
#' combinations of clusters to compare.
#'
#' @param nge The log-normalized gene expression matrix.
#' @param cl The factor with cluster assignments per cell (column of nge).
#' @param deMes The output from \code{\link{fx_calcEScombn}}.
#'
#' @return Differential testing results from Wilcoxon rank sum tests comparing a
#'   gene in each cluster to that gene in every other cluster in a series of
#'   tests. The results are stored as a nested list of dataframes. Each list
#'   element contains a named list of clusters (cluster A). Each of those lists
#'   contains a named list of all the other clusters (cluster B). Each of those
#'   list elements contains a dataframe of four variables, where each sample is
#'   a gene. \code{dDR} is the difference in detection rate of that gene between
#'   the two clusters (DR[A] - DR[B]). \code{logGER} is the log gene expression
#'   ratio calculated by taking the difference in mean expression of the gene
#'   (see \code{\link{meanLogX}} for mean calculation) between the two clusters
#'   (MGE[A] - MGE[B]). \code{Wstat} and \code{pVal} are the test statistic and
#'   the p-value of the Wilcoxon rank sum test. \code{FDR} is the false
#'   discovery rate-corrected p-value of the test.
#'   

fx_calcDEcombn_BP <- function(nge,cl,deMes) {
  combosL <- strsplit(names(deMes),"-")
  message("-- Testing differential expression between clusters --")
  deM_pVal <- BiocParallel::bplapply(seq_along(combosL),function(i)
    apply(nge[rownames(deMes[[i]])[deMes[[i]]$overThreshold],],1,function(X)
      suppressWarnings(unlist(
        wilcox.test(X[cl == combosL[[i]][1]],
                    X[cl == combosL[[i]][2]])[c("statistic","p.value")]
      ))
    )
  )
  for (i in seq_along(deMes)) {
    deMes[[i]][colnames(deM_pVal[[i]]),"Wstat"] <- deM_pVal[[i]]["statistic.W",]
    deMes[[i]][colnames(deM_pVal[[i]]),"pVal"] <- deM_pVal[[i]]["p.value",]
    deMes[[i]][colnames(deM_pVal[[i]]),"FDR"] <- p.adjust(deM_pVal[[i]]["p.value",],"fdr")
  } 
  return(deMes)
}


#' Performs DE testing between pairs of clusters in sCVdata
#'
#' Performs differential gene expression tests between each pairwise combination
#' of cluster in an sCVdata object using the gene expression matrix of input
#' data object. Alternatively, this function can be skipped, and existing DE
#' test results can be assigned directly to the sCVdata object.
#'
#' This function performs Wilcoxon rank sum tests comparing gene expression
#' between the cells of all pairwise combinations of cluster clusters in the
#' input data. Gene expression ratio in log space (\code{logGER}) and
#' differences in detection rates (\code{dDR}) are reported for all genes in the
#' comparison. Genes are tested if they are detected in at least one of the
#' cluster at a higher proportion than \code{Param(sCVd,"DRthresh")}, and both
#' unadjusted p-values and false discovery rates are reported for all genes
#' tested. To help track its progress, this function uses progress bars from
#' \code{pbapply}. To disable these, set
#' \code{\link[pbapply]{pboptions}(type="none")}. To re-enable, set
#' \code{\link[pbapply]{pboptions}(type="timer")}.
#'
#' If using existing DE test results, assign results of differential gene
#' expression tests for all pairwise combinations of clusters in sCVdata to the
#' \code{\link{DEcombn}} slot of the \code{\link{sCVdata}} object. See example
#' and slot documentation.
#'
#' @param sCVd An sCVdata object.
#' @param inD The input dataset. An object of class \code{\link[Seurat]{seurat}}
#'   or \code{\link[SingleCellExperiment]{SingleCellExperiment}}. Other data
#'   classes are not currently supported.
#'   \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#'   for other data objects here!}
#' @param storeAllDE A logical vector of length 1 indicating whether to
#'   calculate and store effect size information for all genes in the comparison
#'   (TRUE), or just those passing the detection rate threshold for the Wilcoxon
#'   rank-sum test (FALSE). Setting this to FALSE will reduce the size of the
#'   output sCVdata object.
#' @param UseBiocParallel Default = FALSE. Very experimental implementation of
#'   BiocParallel for calculations. Not recommended.
#'
#' @return A named list of data frames, one entry for each pairwise combination
#'   of levels in \code{Clusters(sCVd)} (with corresponding name where levels
#'   are separated by '-'). Each entry is data frame containing gene
#'   differential expression stats when comparing the cells of that cluster to
#'   all other cells in the input data. Rows represent genes, and variables
#'   include \code{logGER} (an effect size measure: gene expression ratio in log
#'   space, often referred to as logFC), \code{dDR} (an effect size measure:
#'   difference in detection rate), and \code{FDR} (significance measure: false
#'   discovery rate). Also included are \code{Wstat} and \code{pVal}, the test
#'   statistic and the p-value of the Wilcoxon rank sum test.
#'
#' @seealso \code{\link{CalcSCV}} for wrapper function to calculate all
#'   statistics for an sCVdata object, and \code{\link{fx_calcEScombn}} and
#'   \code{\link{fx_calcDEcombn}} for the internal functions performing the
#'   calculations. Wilcox test is now powered by \code{\link[presto]{wilcoxauc}}
#'   for super speed.
#'
#' @examples
#' \dontrun{
#' ## Example using CalcDEvsRest ##
#' DEcombn(your_sCV_obj) <- CalcDEcombn(sCVd=your_sCV_obj,
#'                                      inD=your_scRNAseq_data_object)
#'
#'
#' ## Example using MAST results from Seurat to replace CalcDEcombn ##
#'
#' MAST_pw <- apply(combn(levels(your_seurat_obj@ident),2),
#'                  MARGIN=2,
#'                  function(X) {
#'                    FindMarkers(your_seurat_obj,
#'                                ident.1=X[1],
#'                                ident.2=X[2],
#'                                logfc.threshold=0,
#'                                min.pct=0.1,
#'                                test.use="MAST",
#'                                latent.vars="nUMI")
#'                  })
# ^ Test DE between every pairwise combination of clusters
#' names(MAST_pw) <- apply(combn(levels(your_seurat_obj@ident),2),2,
#'                         function(X) paste(X,collapse="-"))
#' # ^ Names must be in "X-Y" format
#' 
#' for (i in names(MAST_pw)) {
#'   MAST_pw[[i]]$dDR <- MAST_pw[[i]]$pct.1 - MAST_pw[[i]]$pct.2
#'   # ^ Diff in detect rate (dDR) must be a variable in each dataframe
#'   names(MAST_pw[[i]])[names(MAST_pw[[i]]) == "avg_logFC"] <- "logGER"
#'   # ^ Effect size variable must be named 'logGER'
#'   names(MAST_pw[[i]])[names(MAST_pw[[i]]) == "p_val_adj"] <- "FDR"
#'   # ^ Significance variable must be named 'FDR'
#'   # Note: rownames of each dataframe must be gene names, 
#'   # but FindMarkers should already do this.
#' }
#' DEcombn(your_sCV_obj) <- MAST_pw
#' # ^ Slot MAST results into sCVdata object
#'
#' }
#'
#' @name CalcDEcombn
#'
#' @export
#' 

setGeneric("CalcDEcombn",function(sCVd,inD,storeAllDE,UseBiocParallel)
  standardGeneric("CalcDEcombn"))


#' @describeIn CalcDEcombn Calculate DE between cluster pairs
#' @export

setMethod("CalcDEcombn","sCVdata",
          function(sCVd,inD,storeAllDE=TRUE,UseBiocParallel=FALSE) {
            if (!is(inD)[1] %in% findMethodSignatures(getExpr)) {
              stop(paste("The input data object must be one of:",
                         paste(findMethodSignatures(getExpr),collapse=", "),
                         sep="\n  "))
            }
            if (length(levels(Clusters(sCVd))) <= 1) {
              stop("scClustViz can't calculate differential expression when there's only one cluster.")
            }
            
            deMes <- fx_calcEScombn(cl=Clusters(sCVd),
                                    CGS=ClustGeneStats(sCVd),
                                    DRthresh=Param(sCVd,"DRthresh"))
            if (!storeAllDE) { 
              deMes <- sapply(deMes,function(X) X[X$overThreshold,],simplify=F) 
            }
            if (UseBiocParallel) {
              deMes <- fx_calcDEcombn_BP(nge=getExpr(inD,Param(sCVd,"assayType")),
                                         cl=Clusters(sCVd),
                                         deMes=deMes)
            } else {
              deMes <- fx_calcDEcombn(nge=getExpr(inD,Param(sCVd,"assayType")),
                                      cl=Clusters(sCVd),
                                      deMes=deMes)
            }
            return(deMes)
          })


# CalcSilhouette ----

#' Internal fx to call \code{\link[cluster]{silhouette}} for sCVdata.
#'
#' Internal function. See \code{\link{CalcSilhouette}}.
#' 
#' @param pca The cell embeddings used to calculate clusters.
#' @param cl The cluster assignments.
#' 
#' @return A silhouette object.
#' 

fx_calcSilhouette <- function(pca,cl) {
  return(cluster::silhouette(as.integer(cl),dist(pca)))
}


#' Calculates silhouette widths for all cells in sCVdata
#'
#' Calls \code{\link[cluster]{silhouette}} to calculate silhouette widths (a
#' metric indicating each cell's contribution to cluster cohesion / separation)
#' for all cells in the sCVdata object, using the cell embedding used in
#' clustering (see \code{Param(sCVd,"DRforClust")}).
#'
#' @param sCVd An sCVdata object.
#' @param inD The input dataset. An object of class \code{\link[Seurat]{seurat}}
#'   or \code{\link[SingleCellExperiment]{SingleCellExperiment}}. Other data
#'   classes are not currently supported.
#'   \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#'   for other data objects here!}
#'
#' @return A silhouette object containing the silhouette widths for all cells in
#'   the sCVdata object.
#'
#' @seealso \code{\link{CalcSCV}} for wrapper function to calculate all
#'   statistics for an sCVdata object, \code{\link{fx_calcSilhouette}} for the
#'   internal function this method points to, and
#'   \code{\link[cluster]{silhouette}} for the function doing the calculations.
#'
#' @examples
#' \dontrun{
#' Silhouette(your_sCV_obj) <- CalcSilhouette(sCVd=your_sCV_obj,
#'                                            inD=your_scRNAseq_data_object)
#' }
#'
#' @name CalcSilhouette
#'
#' @export
#' 

setGeneric("CalcSilhouette",function(sCVd,inD) standardGeneric("CalcSilhouette"))


#' @describeIn CalcSilhouette Calculate silhouette widths for sCVdata
#' @export

setMethod("CalcSilhouette",signature("sCVdata"),function(sCVd,inD) {
  if (require(cluster)) {
    if (length(levels(Clusters(sCVd))) <= 1) {
      stop("Silhouette cannot be calculated with a single cluster.")
    }
    fx_calcSilhouette(pca=getEmb(inD,Param(sCVd,"DRforClust")),
                      cl=Clusters(sCVd))
  } else {
    stop(paste("Silhouette could not be calculated because package 'cluster' is missing.",
               "Try 'install.packages(cluster)', then run 'CalcSilhouette()' again.",sep="\n  "))
  }
})


# DEdist ----

#' Internal function for distance between clusters by number of DE genes.
#'
#' Internal function. See \code{\link{DEdist}}.
#'
#' Counts number of DE genes between each pair of clusters from
#' \code{\link{CalcDEcombn}} output.
#'
#' @param deVS List of pairwise DE results from \code{\link{CalcDEcombn}}
#' @param FDRthresh False discovery rate for counting significance.
#'
#' @return A matrix of distances between clusters for each cluster resolution.
#'   Interpretable by \code{\link[stats]{as.dist}} to generate a \code{dist}
#'   object.
#'   

fx_calcDist_numDE <- function(deVS,FDRthresh) {
  deD <- sapply(deVS,function(X) sum(X$FDR <= FDRthresh,na.rm=T))
  names(deD) <- NULL
  return(sapply(unique(unlist(strsplit(names(deVS),"-"))),function(X)
    sapply(unique(unlist(strsplit(names(deVS),"-"))),function(Y) {
      if (X == Y) {
        return(NA)
      } else {
        return(deD[sapply(strsplit(names(deVS),"-"),
                          function(comp) X %in% comp & Y %in% comp)])
      }
    }))) # returns a matrix of distances b/w clusters (for DEneighb)
}


#' Internal function for distance between clusters by DE score.
#'
#' Internal function. See \code{\link{DEdist}}.
#'
#' Sums absolute DE scores between each pair of clusters from
#' \code{\link{CalcDEcombn}} output. DE scores are \code{-log10(FDR)}. Used to
#' determine nearest neighbouring clusters.
#'
#' @param deVS List of pairwise DE results from \code{\link{CalcDEcombn}}
#'
#' @return A matrix of distances between clusters for each cluster resolution.
#'   Interpretable by \code{\link[stats]{as.dist}} to generate a \code{dist}
#'   object.
#'   

fx_calcDist_scoreDE <- function(deVS) {
  d <- vapply(deVS,function(X) {
    temp <- -log10(X$FDR) 
    temp[is.na(temp)] <- 0
    temp[temp == Inf] <- max(temp[temp < Inf]) + 1
    return(sum(temp^2))
    },FUN.VALUE=numeric(1))^0.5
  cb <- strsplit(names(deVS),"-")
  cl <- unique(unlist(cb))
  tempOut <- matrix(nrow=length(cl),ncol=length(cl),dimnames=list(cl,cl))
  for (i in seq_along(d)) {
    tempOut[cb[[i]][1],cb[[i]][2]] <- d[i]
    tempOut[cb[[i]][2],cb[[i]][1]] <- d[i]
  }
  return(tempOut) 
}

#' Internal function for distance between clusters by gene expression statistic.
#'
#' Internal function. See \code{\link{DEdist}}.
#'
#' Calculates euclidean distance in gene expression space between each pair of
#' clusters, using a summary statistic of gene expression for each cluster.
#' Calculated using \code{\link[stats]{dist}}.
#'
#' @param deVS List of pairwise DE results from \code{\link{CalcDEcombn}}
#' @param metric A summary statistic of gene expression per cluster, from
#'   \code{\link{CalcCGS}}
#'
#' @return A matrix of distances between clusters for each cluster resolution.
#'   Interpretable by \code{\link[stats]{as.dist}} to generate a \code{dist}
#'   object.
#'   

fx_calcDist_geneExpr <- function(CGS,metric) {
  if (!metric %in% names(CGS[[1]])) {
    stop(paste("metric must be one of:",paste(names(CGS[[1]]),collapse=", ")))
  }
  tempDist <- as.vector(dist(t(sapply(CGS,function(X) X[[metric]]))))
  tempNames <- apply(combn(names(CGS),2),2,function(X) paste(X,collapse="-"))
  return(sapply(unique(unlist(strsplit(tempNames,"-"))),function(X)
    sapply(unique(unlist(strsplit(tempNames,"-"))),function(Y) {
      if (X == Y) {
        return(NA)
      } else {
        return(tempDist[sapply(strsplit(tempNames,"-"),
                               function(comp) X %in% comp & Y %in% comp)])
      }
    }))) # returns a matrix of distances b/w clusters (for DEneighb)
}

#' Internal function for distance between clusters from cell embedding
#'
#' Internal function. See \code{\link{DEdist}}.
#'
#' Calculates euclidean distance between centroids of each cluster in the cell
#' embedding used for clustering. Calculated using \code{\link[stats]{dist}}.
#'
#' @param deVS List of pairwise DE results from \code{\link{CalcDEcombn}}
#' @param cellEmb A cell embedding with cells in rows and dimensions in columns.
#'
#' @return A matrix of distances between clusters for each cluster resolution.
#'   Interpretable by \code{\link[stats]{as.dist}} to generate a \code{dist}
#'   object.
#'   

fx_calcDist_cellEmb <- function(cl,cellEmb) {
  temp <- as.matrix(dist(apply(cellEmb,2,function(X) tapply(X,cl,mean))))
  temp[temp == 0] <- NA
  return(temp)
}


#' Calculate inter-cluster distances
#'
#' A variety of methods for calculating inter-cluster distances for sCVdata
#' objects, most of which return matrices which can be converted to
#' \code{\link[stats]{dist}} objects.
#'
#' @param sCVd An sCVdata object.
#' @param y One of the following, depending on the desired distance metric:
#'   \describe{
#'     \item{By number of differentially expressed genes:}{A numeric vector of 
#'       length one representing the threshold of false discovery rate 
#'       controlled for when counting significantly differentially expressed 
#'       genes.}
#'     \item{By differential expression significance score:}{Nothing. Only pass
#'       an argument to \code{sCVd}.}
#'     \item{By gene expression summary statistic:}{A character vector of length 
#'       one referring to a column name in \code{ClustGeneStats(sCVd)}: 
#'       \code{MGE}, \code{DR}, or \code{MDGE}.}
#'     \item{By distance between centroids in a cell embedding:}{A matrix 
#'       representing the cell embedding, where cells are rows and dimensions
#'       are columns. See \code{\link{getEmb}} to easily extract this from your
#'       input data object.}
#'   }
#'   
#' @return A matrix of distances between clusters for each cluster resolution.
#'   Interpretable by \code{\link[stats]{as.dist}} to generate a \code{dist}
#'   object.
#'
#' @seealso \code{\link{DEdistNN}} for a method to find nearest neighbouring
#'   clusters from \code{DEdist} results.
#'
#' @examples
#' \dontrun{
#' numDE_bw_clusts <- DEdist(sCVd=your_sCV_obj,0.01)
#' DEdistNN(numDE_bw_clusts)
#' }
#'
#' @name DEdist
#'
#' @export
#' 

setGeneric("DEdist",function(sCVd,y) standardGeneric("DEdist"))


#' @describeIn DEdist By number of differentially-expressed genes (see
#'   \code{\link{fx_calcDist_numDE}})
#' @export

setMethod("DEdist",signature("sCVdata","numeric"),
          function(sCVd,y) fx_calcDist_numDE(deVS=DEcombn(sCVd),FDRthresh=y))


#' @describeIn DEdist By differential expression significance score (see
#'   \code{\link{fx_calcDist_scoreDE}})
#' @export

setMethod("DEdist",signature("sCVdata","missing"),
          function(sCVd,y) fx_calcDist_scoreDE(DEcombn(sCVd)))


#' @describeIn DEdist By gene expression summary statistic (see
#'   \code{\link{fx_calcDist_geneExpr}})
#' @export

setMethod("DEdist",signature("sCVdata","character"),
          function(sCVd,y) fx_calcDist_geneExpr(CGS=ClustGeneStats(sCVd),metric=y))


#' @describeIn DEdist By distance between centroids in a cell embedding (see
#'   \code{\link{fx_calcDist_cellEmb}})
#' @export

setMethod("DEdist",signature("sCVdata","matrix"),
          function(sCVd,y) fx_calcDist_cellEmb(Clusters(sCVd),y))


# DEdistNN ----

#' Determine nearest neighbouring clusters
#'
#' Finds nearest neighbouring clusters from output of \code{\link{DEdist}}.
#'
#' @param d Output of \code{\link{DEdist}}
#'
#' @return A named integer vector where each name is a cluster and element is
#'   the position of its nearest neighbouring cluster in \code{Clusters(sCVd)}.
#'
#' @examples
#' \dontrun{
#' numDE_bw_clusts <- DEdist(sCVd=your_sCV_obj,0.01)
#' DEdistNN(numDE_bw_clusts)
#' }
#'
#' @name DEdistNN
#'
#' @export
#' 

setGeneric("DEdistNN",function(x) standardGeneric("DEdistNN"))


#' @describeIn DEdistNN Nearest neighbours from DEdist vector.
#' @export

setMethod("DEdistNN",signature("numeric"),
          function(x) 
            sapply(unique(unlist(strsplit(names(x),"-"))),
                   function(Y) 
                     as.integer(sub(paste0("-?",Y,"-?"),"",
                                    names(which.min(x[grep(Y,names(x))]))))))


#' @describeIn DEdistNN Nearest neighbours from DEdist matrix.
#' @export

setMethod("DEdistNN",signature("matrix"),
          function(x) apply(x,1,function(X) names(which.min(X))))


# DEneighb ----

#' Internal function to find DE genes between nearest neighbouring clusters
#'
#' Internal function. See \code{\link{DEneighb}}.
#'
#' Identifies the significantly differentially expressed genes between nearest
#' neighbouring clusters, using the pairwise differential gene expression
#' testing results.
#'
#' @param deVS List of pairwise DE results from \code{\link{CalcDEcombn}}.
#' @param NN A vector of nearest neighbours from \code{\link{DEdistNN}}.
#' @param FDRthresh False discovery rate for counting significance.
#'
#' @return A named list of data frames, one entry for each cluster. Each data
#'   frame is the results of the gene expression test comparing the cluster to
#'   its nearest neighbour. See \code{\link{CalcDEcombn}} and
#'   \code{\link{DEcombn}} for details.
#'   


fx_calcNeighb <- function(deVS,NN,FDRthresh) {
  nb <- rbind(paste(names(NN),NN,sep="-"),
              paste(NN,names(NN),sep="-"))
  colnames(nb) <- names(NN)
  nbd <- apply(nb,2,function(X) X %in% names(deVS))
  nb <- nb[nbd]
  nbd <- apply(nbd,2,which)
  nbd[nbd == 2] <- -1
  if (missing(FDRthresh)) { FDRthresh <- 1 }
  deN <- sapply(seq_along(nb),function(i) {
    temp <- which(deVS[[nb[i]]]$FDR <= FDRthresh &
                    deVS[[nb[i]]]$logGER * nbd[i] > 0)
    out <- deVS[[nb[i]]][temp,names(deVS[[nb[i]]]) != "overThreshold"]
    names(out) <- paste(names(out),paste(names(NN),NN,sep="-")[i],sep="_")
    return(out)
  },simplify=F)
  names(deN) <- names(NN)
  return(deN)
}


#' Find differentially expressed genes between nearest neighbouring clusters
#'
#' Identifies the significantly differentially expressed genes between nearest
#' neighbouring clusters, using the pairwise differential gene expression
#' testing results.
#'
#' @param sCVd An sCVdata object.
#' @param FDRthresh False discovery rate for counting significance.
#'
#' @return A named list of data frames, one entry for each cluster. Each data
#'   frame is the results of the gene expression test comparing the cluster to
#'   its nearest neighbour. See \code{\link{CalcDEcombn}} and
#'   \code{\link{DEcombn}} for details.
#'
#' @name DEneighb
#'
#' @export
#' 

setGeneric("DEneighb",function(sCVd,FDRthresh) standardGeneric("DEneighb"))


#' @describeIn DEneighb Find DE genes between nearest neighbours from sCVdata
#' @export

setMethod("DEneighb","sCVdata",function(sCVd,FDRthresh) fx_calcNeighb(deVS=DEcombn(sCVd),
                                                                      NN=DEdistNN(DEdist(sCVd)),
                                                                      FDRthresh=FDRthresh))


# DEmarker ----

#' Internal function to find marker genes for each cluster
#'
#' Internal function. See \code{\link{DEmarker}}.
#'
#' Identifies the significantly positively differentially expressed genes
#' between a cluster and all other clusters, using the pairwise differential
#' gene expression testing results.
#'
#' @param deVS List of pairwise DE results from \code{\link{CalcDEcombn}}.
#' @param FDRthresh False discovery rate for counting significance.
#'
#' @return A named list of data frames, one entry for each cluster. Each data
#'   frame is the combined results of the set of pairwise gene expression tests
#'   comparing the cluster to each other cluster in the data. See
#'   \code{\link{CalcDEcombn}} and \code{\link{DEcombn}} for details.
#'   

fx_calcMarker <- function(deVS,FDRthresh) {
  combosL <- sapply(unique(unlist(strsplit(names(deVS),"-"))),
                    function(clust) {
                      temp <- sapply(strsplit(names(deVS),"-"),
                                     function(comp) 
                                       which(comp == clust))
                      out <- unlist(temp)
                      out[out == 2] <- -1
                      names(out) <- which(sapply(temp,length) > 0)
                      return(out)
                    },
                    simplify=F)
  if (missing(FDRthresh)) { FDRthresh <- 1 }
  mNames <- sapply(combosL,function(comp) {
    cn <- as.integer(names(comp))
    if ("Wstat" %in% colnames(deVS[[1]])) {
      return(Reduce(intersect,sapply(seq_along(comp),function(i) {
        tempW <- deVS[[cn[i]]]$Wstat - deVS[[cn[i]]]$Wstat[which.max(deVS[[cn[i]]]$pVal)]
        rownames(deVS[[cn[i]]])[which(deVS[[cn[i]]]$FDR <= FDRthresh & tempW * comp[i] > 0)]
      },simplify=F)))
    } else {
      return(Reduce(intersect,sapply(seq_along(comp),function(i) {
        rownames(deVS[[cn[i]]])[which(deVS[[cn[i]]]$FDR <= FDRthresh & 
                                             deVS[[cn[i]]]$logGER * comp[i] > 0)]
      },simplify=F)))
    }
  },simplify=F)
  deM <- sapply(seq_along(mNames),function(i) {
    do.call(cbind,lapply(names(combosL[[i]]),function(l) {
      L <- as.integer(l)
      temp <- deVS[[L]][mNames[[i]],names(deVS[[L]]) != "overThreshold"]
      if (combosL[[i]][l] == -1) {
        temp$logGER <- temp$logGER * combosL[[i]][l]
        temp$dDR <- temp$dDR * combosL[[i]][l]
        tempN <- paste(strsplit(names(deVS)[L],split="-")[[1]][c(2,1)],collapse="-")
        names(temp) <- paste(names(temp),tempN,sep="_")
      } else {
        names(temp) <- paste(names(temp),names(deVS)[L],sep="_")
      }
      return(temp)
    }))
  },simplify=F)
  names(deM) <- names(mNames)
  return(deM)
}


#' Find marker genes for each cluster
#'
#' Identifies the significantly positively differentially expressed genes
#' between a cluster and all other clusters, using the pairwise differential
#' gene expression testing results.
#'
#' @param sCVd An sCVdata object.
#' @param FDRthresh False discovery rate for counting significance.
#'
#' @return A named list of data frames, one entry for each cluster. Each data
#'   frame is the combined results of the set of pairwise gene expression tests
#'   comparing the cluster to each other cluster in the data. See
#'   \code{\link{CalcDEcombn}} and \code{\link{DEcombn}} for details.
#'
#' @name DEmarker
#'
#' @export
#' 

setGeneric("DEmarker",function(sCVd,FDRthresh) standardGeneric("DEmarker"))


#' @describeIn DEmarker Find marker genes from sCVdata
#' @export

setMethod("DEmarker","sCVdata",function(sCVd,FDRthresh) fx_calcMarker(deVS=DEcombn(sCVd),
                                                                      FDRthresh=FDRthresh))

