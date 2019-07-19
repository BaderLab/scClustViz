#' Computes the mean of log-scaled values
#' 
#' Computes the arithmetic mean in linear space of log-scaled values, and returns the result
#' in the same log scale.
#' 
#' Generally a pseudocount of 1 is added to log-scaled values to prevent +/-Inf results. 
#' However, adding a pseudocount of 1 to the log-scaled mean prior to gene expression ratio 
#' calculations in log space skews the result quite dramatically, so instead we add a small 
#' pseudocount to avoid +/- inf results when means are zero, without the same skewing. Adding 
#' a very small (ie 1e-99) number means that means of zero get set to a large negative log-mean, 
#' when it might be more appropriate to have those values fall closer to the smallest non-zero 
#' log-mean. By using a pseudocount of 1 / number of samples in the experiment, we ensure that 
#' log(zero) is smaller than any non-zero log-mean, while still being in the same ballpark.
#' 
#' @param data A numeric vector in log scale for which the arithmetic mean in linear scale 
#' is to be calculated and returned in the same log scale. 
#' @param ncell The number of samples (cells assuming an scRNAseq experiment) in the data. 
#' This is not necessarily the same as the length of \code{data}. Since its inverse will 
#' be used as the pseudocount value for all calculations, it should be consistent for all 
#' calls of \code{mean.logX}.
#' @param ex The log base of the data. Traditionally gene expression data is represented 
#' in base 2, although some methods (ie. Seurat's normalization scheme) use the natural log. 
#' (default=2)
#' @param pc The pseudocount used when converting the data to log-scale initially. (default=1)
#' 
#' @export
#' 

meanLogX <- function(data,ncell,ex=2,pc=1) { 
  log(mean(ex^data - pc) + 1/ncell,base=ex)
}


#' Computes the cosine similarity between two vectors
#'
#' Computes the cosine similarity between two vectors.
#' 
#' @references \url{https://stats.stackexchange.com/q/31573}/
#'
#' @param A One vector in the pair to compare.
#' @param B The other vector in the pair to compare.
#'
#' @return A numeric value, the cosine similarity between vectors \code{A} and
#'   \code{B}
#' 

cosineSim <- function(A,B) sum(A*B)/sqrt(sum(A^2)*sum(B^2))


#' Automatically determines keytype for AnnotationDb lookup
#'
#' Compares rownames of gene expression matrix with keys for all keytypes in an
#' AnnotationDb object, and returns the keytype with the best match.
#'
#' If the rownames for the gene expression matrix are not official gene symbols
#' (determined by less than 80 percent of rownames matching keys in the SYMBOL keytype
#' of the AnnotationDb object), this function compares rownames of gene
#' expression matrix with keys for all keytypes in an AnnotationDb object, and
#' returns the keytype with the best match.
#'
#' @param nge The gene expression matrix, i.e. \code{getExpr(yourDataObject)}.
#' @param annotationDB The AnnotationDb object.
#'
#' @return A character value, the \code{keytype} of the AnnotationDb object, to
#'   be passed to \code{map2symbol}.
#' 
#' @seealso \code{\link{map2symbol}}
#' 
#' @export

findKeyType <- function(nge,annotationDB) {
  rownameKeytype <- "SYMBOL"
  if (sum(rownames(nge) %in% keys(annotationDB,rownameKeytype)) / nrow(nge) < 0.8) {
    warning(paste("Less than 80% of rownames map to official gene symbols.",
                  "Automatically determining keytype from rownames..."))
    temp_keyMatch <- pbapply::pbsapply(AnnotationDbi::keytypes(annotationDB),function(X) 
      sum(rownames(nge) %in% AnnotationDbi::keys(annotationDB,X)))
    rownameKeytype <- names(which.max(temp_keyMatch))
    warning(paste0("Keytype '",rownameKeytype,"' matched ",
                 max(temp_keyMatch),"/",nrow(nge)," rownames."))
  }
  return(rownameKeytype)
}


#' Maps gene identifiers to official gene symbol
#'
#' Maps gene identifiers in the gene expression matrix to official gene symbols
#' for ease of data interpretation in scClustViz UI.
#'
#' If rownames are already official gene symbols, returns null.  In this case,
#' \code{addCellMarkersToCGS} and \code{labelCellTypes} will proceed under the
#' assumption that the gene identifiers in the input data are already official
#' gene symbols.
#'
#' @param nge The gene expression matrix, i.e. \code{getExpr(yourDataObject)}.
#' @param annotationDB The AnnotationDb object.
#' @param rownameKeytype The keytype of the gene identifiers. Can be determined
#'   automatically using \code{findKeyType}.
#'
#' @return A named character vector, the output of \code{AnnotationDbi::mapIds}
#'   where values are gene symbols and names are the input gene identifiers.
#'
#' @seealso \code{\link{findKeyType}},\code{\link[AnnotationDbi]{mapIds}},
#'   \code{\link{addCellMarkersToCGS}}, and \code{\link{labelCellTypes}}.
#' 
#' @export

map2symbol <- function(nge,annotationDB,rownameKeytype) {
  if (rownameKeytype != "SYMBOL") {
    return(AnnotationDbi::mapIds(x=annotationDB,
                                 keys=rownames(nge),
                                 column="SYMBOL",
                                 keytype=rownameKeytype,
                                 multiVals="first"))
  } else {
    return(NULL)
  }
}


#' Add gene symbols and cell type marker information to cluster gene statistics
#'
#' Adds four new variables to the cluster-wise gene statistics dataframe of the
#' scClustViz data object. Official gene symbols are added as variable
#' \code{genes}. The remaining variables are used in
#' \code{\link{plot_clusterGenes_markers}} to plot cell type marker genes.
#' Variables \code{cMu} and \code{cMs} are logical vectors indicating genes that
#' are unique to and shared across cell type markers respectively. Variable
#' \code{overCut} indicates which genes should be labelled in the plot.
#'
#' @param sCV An object of class \code{sCVdata}.
#' @param cellMarkersU Derived from the \code{cellMarkers} argument to
#'   \code{\link{runShiny}}. A list of the unique gene symbols for each cell
#'   type in \code{cellMarkers}.
#' @param cellMarkersS Derived from the \code{cellMarkers} argument to
#'   \code{\link{runShiny}}. A list of the gene symbols common to two or more
#'   cell types in \code{cellMarkers}. Each entry is named for the indicies of
#'   \code{cellMarkers} that share the gene.
#' @param symbolMap The output of \code{\link{map2symbol}}.
#'
#' @return The \code{sCVdata} object with the new variables added to
#'   \code{ClustGeneStats}.
#'
#' @seealso \code{\link{findKeyType}}, \code{\link{map2symbol}}, and
#'   \code{\link[AnnotationDbi]{mapIds}}.
#'   
#' @export

addCellMarkersToCGS <- function(sCV,cellMarkersU,cellMarkersS,symbolMap) {
  if (is.null(ClustGeneStats(sCV))) {
    stop("ClustGeneStats(sCV) cannot be NULL. Run CalcAllSCV or CalcSCV before calling runShiny.")
  }
  ClustGeneStats(sCV) <- sapply(ClustGeneStats(sCV),function(CGS) {
    if (!is.null(symbolMap)) {
      CGS$genes <- symbolMap[rownames(CGS)]
      CGS$genes[is.na(CGS$genes)] <- rownames(CGS)[is.na(CGS$genes)]
    } else {
      CGS$genes <- rownames(CGS)
    }
    CGS$cMu <- 
      rownames(CGS) %in% unlist(cellMarkersU) |
      CGS$genes %in% unlist(cellMarkersU)
    CGS$cMs <- (rownames(CGS) %in% unlist(cellMarkersS) | 
                  CGS$genes %in% unlist(cellMarkersS))
    CGS$overCut <- CGS$MGE > mean(CGS$MGE)
    return(CGS)
  },simplify=F)
  return(sCV)
}


#' scClustViz helper fx: Add predicted cell type names to cluster labels
#'
#' A bare-bones method of predicting cell types from marker genes.
#'
#' Assigns cell type labels to each cluster based on the rank of median gene
#' expression of marker genes for each cell type. There are many more
#' intelligent methods for cell type prediction out there.  See
#' \href{https://doi.org/10.1101/521914}{CellAssign}, for example.
#'
#' @param sCV An object of class \code{sCVdata}.
#' @param cellMarkers The \code{cellMarkers} argument from
#'   \code{\link{runShiny}}. A list of marker genes for expected cell types.
#' @param symbolMap Default=NULL. The output of \code{\link{map2symbol}}. If the
#'   rownames (gene identifiers) of your input data object match the gene
#'   identifiers used in your \code{cellMarkers} list, you can leave this as
#'   \code{NULL}, since no gene identifier mapping needs to be performed.
#'
#' @return Returns the sCVdata object with an added attribute
#'   '\code{ClusterNames}' to \code{Clusters(sCV)} containing the assigned cell
#'   type names for each cluster. Stores the \code{cellMarkers} list as an
#'   attribute in \code{Clusters(sCV)}. Also adds four new variables to
#'   \code{ClustGeneStats(sCV)}: Official gene symbols are added as variable
#'   \code{genes}. The remaining variables are used in
#'   \code{\link{plot_clusterGenes_markers}} to plot cell type marker genes in
#'   the Shiny app (see \code{\link{runShiny}}). Variables \code{cMu} and
#'   \code{cMs} are logical vectors indicating genes that are unique to and
#'   shared across cell type markers respectively. Variable \code{overCut}
#'   indicates which genes should be labelled in the plot.
#'
#' @export

labelCellTypes <- function(sCV,cellMarkers,symbolMap=NULL) {
  if (missing(cellMarkers)) {
    cellMarkers <- list() 
  }
  if (!is.list(cellMarkers)) {
    stop("cellMarkers must be a list where each entry is named for a cell type",
         "and is a character vector of gene names for cell type markers.")
  }
  if (length(cellMarkers) < 1) {
    cellMarkersS <- cellMarkersU <- list()
  } else {
    cellMarkersS <- apply(combn(seq_along(cellMarkers),2),2,
                          function(X) do.call(intersect,unname(cellMarkers[X])))
    try(names(cellMarkersS) <- apply(combn(seq_along(cellMarkers),2),2,
                                     function(X) paste(X,collapse="&")),silent=T)
    cellMarkersS <- cellMarkersS[sapply(cellMarkersS,length) > 0]
    cellMarkersU <- lapply(cellMarkers,function(X) X[!X %in% unlist(cellMarkersS)])
  }
  
  sCV <- addCellMarkersToCGS(sCV=sCV,
                             cellMarkersU=cellMarkersU,
                             cellMarkersS=cellMarkersS,
                             symbolMap=symbolMap)
  if (!is.null(attr(Clusters(sCV),"ClusterNames"))) {
    return(sCV)
  } 
  
  if (length(cellMarkers) < 1) {
    attr(Clusters(sCV),"ClusterNames") <- vapply(ClustGeneStats(sCV),
                                                 FUN.VALUE=character(1),
                                                 function(X) return(""))
  } else if (
    if (!is.null(symbolMap)) {
      !any(unlist(cellMarkers) %in% c(symbolMap,names(symbolMap)))
    } else {
      !any(unlist(cellMarkers) %in% rownames(ClustGeneStats(sCV)[[1]]))
    }
  ) {
    warning(paste("None of the provided cellMarkers are found in the data",
                  "(check your gene IDs against rownames in your data)."))
    attr(Clusters(sCV),"ClusterNames") <- vapply(ClustGeneStats(sCV),
                                FUN.VALUE=character(1),
                                function(X) return(""))
  } else {
    temp <- vapply(ClustGeneStats(sCV),
                   FUN.VALUE=character(1),
                   function(Z) 
                     names(which.max(sapply(cellMarkers,function(X) 
                       median(Z$MGE[Z$genes %in% X])))))
    temp[names(temp) == "Unselected"] <- "Unselected"
    attr(Clusters(sCV),"ClusterNames") <- temp
    attr(Clusters(sCV),"cellMarkers") <- cellMarkers
  }
  return(sCV)
}
