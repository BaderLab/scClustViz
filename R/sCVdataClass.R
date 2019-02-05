# sCVparams ----

#' Parameters used in the sCVdata class and methods.
#'
#' An S4 class storing parameters used in the sCVdata class and methods. Access
#' this object from the containing \code{\link{sCVdata}} object with
#' \code{\link{Param}}. Slots are accessed using
#' \code{\link{Param}(sCVdata,slotName)}
#'
#' @slot exponent A length-one numeric vector representing the base of the
#'   log-normalized gene expression data to be processed. Generally gene
#'   expression data is transformed into log2 space when normalizing (set this
#'   to 2), though \code{Seurat} uses the natural log (set this to exp(1)).
#' @slot pseudocount A length-one numeric vector representing the pseudocount
#'   added to the data during log-normalization to avoid log(0) errors.
#'   Generally this is 1.
#' @slot DRthresh A length-one numeric vector between 0 and 1 representing the
#'   detection rate threshold for inclusion of a gene in the differential
#'   expression testing. A gene will be included if it is detected in at least
#'   this proportion of cells in at least one of the clusters being compared.
#'   Commonly set to 0.1.
#' @slot DRforClust A length-one character vector representing the
#'   dimensionality reduction method used as the input for clustering. This is
#'   commonly PCA, and should correspond to the slot name of the cell embedding
#'   in your input data - either the \code{type} argument in
#'   \code{\link[SingleCellExperiment]{reducedDim}(x,type)} or the
#'   \code{reduction.type} argument in
#'   \code{\link[Seurat]{GetDimReduction}(object,reduction.type)} (v2) or
#'   \code{reduction} in \code{\link[Seurat]{Embeddings}(object,reduction)}.
#'
#' @seealso \code{\link{sCVdata}} for containing class.
#'
#' @export
#' 

sCVparams <- setClass(Class="sCVparams",
                      slots=c(exponent="numeric",
                              pseudocount="numeric",
                              DRthresh="numeric",
                              DRforClust="character"))
setValidity("sCVparams",function(object) {
  if (length(object@exponent) > 1) {
    return(paste("exponent refers to the log base of the normalized data",
                 "(generally 2, Seurat uses exp(1)) and should be a single value."))
  } 
  if (length(object@pseudocount) > 1) {
    return(paste("pseudocount refers to the pseudocount used in the",
                 "log-normalization (generally 1) and should be a single value."))
  } 
  if (length(object@DRthresh) > 1) {
    if (object@DRthresh[1] > 1 | object@DRthresh[1] <= 0 ) {
      return(paste("DRthresh refers to detection rate threshold for inclusion",
                   "of a gene in the differential expression testing. A gene will be",
                   "included if it is detected in at least this proportion of cells in",
                   "at least one of the clusters being compared.",
                   "DRthresh should be a single value between 0 and 1"))
    }
  }
  if (length(object@DRforClust) > 1) {
    return(paste("DRforThresh refers to the dimensionality reduction method used",
                 "prior to clustering, and should be a character vector of length 1",
                 "that would be used as the 'type' argument in",
                 "SingleCellExperiment::reducedDim(x,type) or the 'reduction.type'",
                 "argument in Seurat::GetDimReduction(object,reduction.type)"))
  }
})


# sCVdata ----

#' Analysis results for the scClustViz app
#'
#' An S4 class to store analysis results for the scClustViz app. sCVdata objects
#' should be generated using the function \code{\link{CalcSCV}}, which at
#' minimum populates the slots \code{Clusters}, \code{ClustGeneStats}, and
#' \code{params}. The remaining slots can be optionally populated - see slot
#' entries for details. For efficiency, multiple sCVdata objects can be
#' generated for multiple cluster solutions using the single command
#' \code{\link{CalcAllSCV}}.
#'
#' @slot Clusters A named factor representing cluster assignments for every
#'   cell. Accessed with \code{\link{Clusters}}. Length should be equal to
#'   number of cells in input data, and names should match colnames of gene
#'   expression matrix (generally cell barcodes). Levels can be cluster numbers
#'   or cluster names.
#' @slot ClustGeneStats A named list of data frames, one entry for each level in
#'   \code{Clusters} (with corresponding name). Accessed with
#'   \code{ClustGeneStats}. Each entry is data frame containing gene summary
#'   statistics for the cluster. Each data frame has the same number of rows as
#'   the input gene expression matrix, where each row represents the results for
#'   that gene (and shares its name).  The three variables in the data frame are
#'   summary statistics. Detection Rate (\code{DR}) refers to the proprotion of
#'   cells in the cluster in which gene cDNA was detected (a gene expression
#'   value > 0). Mean Detected Gene Expression (\code{MDGE}) refers to the mean
#'   normalized gene expression in cells from the cluster in which the gene was
#'   detected (i.e. non-zero mean). Mean Gene Expression (\code{MGE}) is the
#'   mean normalized gene expression of the gene in the cluster.
#' @slot DEvsRest A named list of data frames, one entry for each level in
#'   \code{Clusters} (with corresponding name). Accessed with \code{DEvsRest}.
#'   Each entry is data frame containing gene differential expression stats when
#'   comparing the cells of that cluster to all other cells in the input data,
#'   as calculated by \code{CalcDEvsRest}. Rows represent genes, and variables
#'   \strong{must include} \code{logGER} (an effect size measure: gene
#'   expression ratio in log space, often referred to as logFC) and \code{FDR}
#'   (significance measure: false discovery rate). This slot can be populated
#'   using the function \code{\link{CalcDEvsRest}}, which can either calculate
#'   differential expression for all clusters, or take an appropriately
#'   formatted list of precomputed values. See \code{\link{CalcDEvsRest}}
#'   documentation for details.
#' @slot DEcombn A named list of data frames, one entry for each combination of
#'   pairs of levels in \code{Clusters} (named as ClusterA-ClusterB). Accessed
#'   with \code{DEcombn}. Each entry is a data frame containing gene
#'   differential expression stats when comparing the two clusters, as
#'   calculated by \code{\link{CalcDEvsCombn}}. Rows represent genes, and
#'   variables \strong{must include} \code{logGER} (an effect size measure: gene
#'   expression ratio in log space, often referred to as logFC) and \code{FDR}
#'   (significance measure: false discovery rate). This slot can be populated
#'   using the function \code{\link{CalcDEvsCombn}}, which can either calculate
#'   differential expression for all combinations of clusters, or take an
#'   appropriately formatted list of precomputed values. See
#'   \code{\link{CalcDEvsCombn}} documentation for details.
#' @slot Silhouette An object of class Silhouette defining a cluster
#'   cohesion/separation statistic (silhouette width) for every cell. Accessed
#'   with \code{Silhouette}. This slot can be populated using the function
#'   \code{\link{CalcSilhouette}}, which is a wrapper to
#'   \code{\link[cluster]{silhouette}} with distance calcuated using the same
#'   reduced dimensional cell embedding as was used for clustering (see
#'   \code{Param(sCVdata,"DRforClust")}).
#' @slot params An object of class \code{\link{sCVparams}} containing the
#'   parameters relevant to the data in this object. Accessed with \code{Param}.
#'   Slots are accessed using \code{\link{Param}(sCVdata,slotName)}.
#'
#' @export
#' 

sCVdata <- setClass(Class="sCVdata",
                    slots=c(Clusters="factor", #need to strip '-' from cluster names
                            ClustGeneStats="list", #list of length(getClusters(sCVdata)) containing data frames 
                            DEvsRest="list", #list of length(getClusters(sCVdata)) containing data frames
                            DEcombn="list", #list of length(getClusters(sCVdata)) containing lists of length(getClusters(sCVdata))-1 containing data frames 
                            Silhouette="ANY",
                            params="sCVparams"))

setValidity("sCVdata",function(object) {
  if (length(object@Clusters) > 0) {
    if (is.null(names(object@Clusters))) {
      return("Clusters should be a named factor where names are cells (colnames) of the input data.")
    }
    if (any(grepl("-",levels(object@Clusters)))) {
      return("Cluster names cannot contain '-'.")
    }
  }
  if (length(object@ClustGeneStats) > 0) {
    if (!identical(names(object@ClustGeneStats),levels(object@Clusters)) |
        !all(sapply(object@ClustGeneStats,is.data.frame))) {
      return(paste("ClustGeneStats should be built using function 'calcClustGeneStats'.",
            "Expected format is a named list with a data frame for each level in @Clusters."))
    }
  }
  if (length(object@DEvsRest) > 0) {
    if (!identical(names(object@DEvsRest),levels(object@Clusters)) |
        !all(sapply(object@DEvsRest,is.data.frame))) {
      return(paste("DEvsRest should be a named list with data frames for each level in @Clusters."))
    }
    if (!all(sapply(object@DEvsRest,
                    function(X) 
                      all(c("logGER","FDR") %in% names(X))))) {
      return("All DEvsRest data frames must contain variables 'logGER' and 'FDR'.")
    }
  }

  if (length(object@DEcombn) > 0) {
    if (length(object@DEcombn) != choose(length(levels(object@Clusters)),2) |
        !all(sapply(object@DEcombn,is.data.frame))) {
      return(paste("DEcombn should be a named list with data frames",
                   "for pairwise combinations of levels in @Clusters."))
    }
    tempNames <- apply(combn(levels(object@Clusters),2),2,
                       function(X) c(paste(X,collapse="-"),
                                     paste(X[2:1],collapse="-")))
    if (!all(sapply(seq_along(object@DEcombn),function(X) 
      names(object@DEcombn)[X] %in% tempNames[,X]))) {
      return(paste("Each entry name in DEcombn should be the names of the pair",
              "of clusters (levels of @Clusters) separated by '-'.",
              "See ?CalcDEcombn for example.",sep="\n"))
    }
    if (!all(sapply(object@DEcombn,
                    function(X) 
                      all(c("logGER","dDR","FDR") %in% names(X))))) {
      return("All DEcombn data frames must contain variables 'logGER', 'dDR', and 'FDR'.")
    }
  }
  
  if (!is.null(object@Silhouette)) {
    if (is(object@Silhouette)[1] != "silhouette") {
      return("Silhouette should be of the class 'silhouette' as returned by cluster::silhouette.")
    }
  }
  validObject(object@params)
})


# ^ slot getters ----------

setGeneric("Clusters",function(sCVd) standardGeneric("Clusters"))
#' @describeIn sCVdata Access Clusters slot
#' @aliases Clusters
#' @export
setMethod("Clusters","sCVdata",function(sCVd) sCVd@Clusters)
setGeneric("Clusters<-",function(sCVd,value) standardGeneric("Clusters<-"))
#' @describeIn sCVdata Assign Clusters slot
#' @export
setReplaceMethod("Clusters","sCVdata",
                 function(sCVd,value) initialize(sCVd,Clusters=value))

setGeneric("ClustGeneStats",function(sCVd) standardGeneric("ClustGeneStats"))
#' @describeIn sCVdata Access ClustGeneStats slot
#' @aliases ClustGeneStats
#' @export
setMethod("ClustGeneStats","sCVdata",function(sCVd) sCVd@ClustGeneStats)
setGeneric("ClustGeneStats<-",function(sCVd,value) standardGeneric("ClustGeneStats<-"))
#' @describeIn sCVdata Assign ClustGeneStats slot
#' @export
setReplaceMethod("ClustGeneStats","sCVdata",
                 function(sCVd,value) initialize(sCVd,ClustGeneStats=value))

setGeneric("DEvsRest",function(sCVd) standardGeneric("DEvsRest"))
#' @describeIn sCVdata Access DEvsRest slot
#' @aliases DEvsRest
#' @export
setMethod("DEvsRest","sCVdata",function(sCVd) sCVd@DEvsRest)
setGeneric("DEvsRest<-",function(sCVd,value) standardGeneric("DEvsRest<-"))
#' @describeIn sCVdata Assign DEvsRest slot - see slot details
#' @export
setReplaceMethod("DEvsRest","sCVdata",
                 function(sCVd,value) initialize(sCVd,DEvsRest=value))

setGeneric("DEcombn",function(sCVd) standardGeneric("DEcombn"))
#' @describeIn sCVdata Access DEcombn slot
#' @aliases DEcombn
#' @export
setMethod("DEcombn","sCVdata",function(sCVd) sCVd@DEcombn)
setGeneric("DEcombn<-",function(sCVd,value) standardGeneric("DEcombn<-"))
#' @describeIn sCVdata Assign DEcombn slot - see slot details
#' @export
setReplaceMethod("DEcombn","sCVdata",
                 function(sCVd,value) initialize(sCVd,DEcombn=value))

setGeneric("Silhouette",function(sCVd) standardGeneric("Silhouette"))
#' @describeIn sCVdata Access Silhouette slot
#' @aliases Silhouette
#' @export
setMethod("Silhouette","sCVdata",function(sCVd) sCVd@Silhouette)
setGeneric("Silhouette<-",function(sCVd,value) standardGeneric("Silhouette<-"))
#' @describeIn sCVdata Assign Silhouette slot
#' @export
setReplaceMethod("Silhouette","sCVdata",
                 function(sCVd,value) initialize(sCVd,Silhouette=value))

setGeneric("Param",function(sCVd,param) standardGeneric("Param"))
#' @describeIn sCVdata Access Param slot (see \code{\link{sCVparams}})
#' @aliases Param
#' @export
setMethod("Param","sCVdata",
          function(sCVd,param) {
            if (missing("param")) {
              slot(sCVd,"params")
            } else {
              if (param %in% slotNames(slot(sCVd,"params"))) {
                slot(sCVd@params,param)
              } else {
                stop(paste("param must be one of:",
                           paste(slotNames(slot(sCVd,"params")),
                                 collapse=", ")))
              }
            }
          })
