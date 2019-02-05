# Generics & methods for loading data from various single-cell data classes.

#' @include sCVdataClass.R
NULL


# Generics ----

# ^ getExpr ----

#' Get gene expression matrix from input data object
#'
#' Extract the gene expression matrix from a single-cell data object containing
#' the input data for scClustViz visualization.
#'
#' This is a wrapper function to the relevant class's normalized data slot
#' accessor method. Currently supported input object classes:
#' \itemize{
#'   \item Class \code{\link[Seurat]{seurat}} accessed by 
#'     \code{\link[Seurat]{GetAssayData}(x)}.
#'   \item Class \code{\link[SingleCellExperiment]{SingleCellExperiment}} 
#'     accessed by \code{\link[SingleCellExperiment]{logcounts}(x)}.
#' }
#' \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#' for other data objects here!}
#' 
#' @param x The single-cell data object.
#' @name getExpr
#' @export
#' 
setGeneric("getExpr",function(x) standardGeneric("getExpr"))


# ^ getMD ----

#' Get metadata from input data object
#'
#' Extract the cell metadata data frame from a single-cell data object
#' containing the input data for scClustViz visualization.
#'
#' This is a wrapper function to the relevant class's cell metadata slot
#' accessor / assignment method. Currently supported input object classes:
#' \itemize{
#'   \item Class \code{\link[Seurat]{seurat}} accessed by 
#'     \code{slot(x,"meta.data")}.
#'   \item Class \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'     accessed by \code{\link[SingleCellExperiment]{colData}(x)}.
#' }
#' \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#' for other data objects here!}
#'
#' @param x The single-cell data object.
#' @name getMD
#' @export
#' 
setGeneric("getMD",function(x) standardGeneric("getMD"))

#' @rdname getMD
#' @export
setGeneric("getMD<-",function(x,value) standardGeneric("getMD<-"))


# ^ getEmb ----

#' Get cell embeddings from input data object
#'
#' Extract the cell embedding coordinates from the dimensionality reduction
#' results stored in the single-cell data object containing the input data for
#' scClustViz visualization.
#'
#' This is a wrapper function to the relevant class's cell metadata slot
#' accessor / assignment method. Currently supported input object classes:
#' \itemize{
#'   \item Class \code{\link[Seurat]{seurat}} accessed by 
#'     \code{\link[Seurat]{GetCellEmbeddings}(x,DRtype)}.
#'   \item Class \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'     accessed by \code{\link[SingleCellExperiment]{reducedDim}(x,DRtype)}.
#' }
#' \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#' for other data objects here!}
#'
#' @param x The single-cell data object.
#' @param DRtype The slot name (generally acronym of the DR method) where the
#'   cell coordinates are stored, case-insensitive.
#' @name getEmb
#' @export
#' 
setGeneric("getEmb",function(x,DRtype) standardGeneric("getEmb"))


# Methods ----

# ^ Seurat ----
suppressMessages(
  setMethod("getExpr","seurat",
            function(x) Seurat::GetAssayData(x))
)

suppressMessages(
  setMethod("getMD","seurat",
            function(x) slot(x,"meta.data"))
)

suppressMessages(
  setReplaceMethod("getMD","seurat",
                   function(x,value) initialize(x,meta.data=value))
)

suppressMessages(
  setMethod("getEmb","seurat",function(x,DRtype) {
    if (tolower(DRtype) %in% names(slot(x,"dr"))) {
      Seurat::GetCellEmbeddings(x,tolower(DRtype))
    } else {
      stop(paste(paste0("DRtype '",DRtype,"' not found."),
                 "The following cell embeddings are available in this object:",
                 paste0(names(slot(x,"dr")),collapse=", "),sep="\n  "))
    }
  })
)

# ^ SingleCellExperiment ----
suppressMessages(
  setMethod("getExpr","SingleCellExperiment",
            function(x) SingleCellExperiment::logcounts(x))
)

suppressMessages(
  setMethod("getMD","SingleCellExperiment",
            function(x) SingleCellExperiment::colData(x))
)

suppressMessages(
  setReplaceMethod("getMD","SingleCellExperiment",
                   function(x,value) initialize(x,colData=value))
)

suppressMessages(
  setMethod("getEmb","SingleCellExperiment",function(x,DRtype) { 
    if (any(grepl(DRtype,SingleCellExperiment::reducedDimNames(x),ignore.case=T))) {
      SingleCellExperiment::reducedDim(x,grep(DRtype,
                                              SingleCellExperiment::reducedDimNames(x),
                                              ignore.case=T,value=T))
    } else {
      stop(paste(paste0("DRtype '",DRtype,"' not found."),
                 "The following cell embeddings are available in this object:",
                 paste0(SingleCellExperiment::reducedDimNames(x),collapse=", "),
                 sep="\n  "))
    }
  })
)