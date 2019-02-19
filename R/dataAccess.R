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
#'   \item Class \code{\link[Seurat]{seurat}/\link[Seurat]{Seurat}} stored in
#'     \code{x@data} or \code{x@assays[[assayType]]@data}, 
#'     depending on Seurat object version.
#'   \item Class \code{\link[SingleCellExperiment]{SingleCellExperiment}} 
#'     accessed by \code{\link[SummarizedExperiment]{assay}(x,assayType)}.
#' }
#' \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#' for other data objects here!}
#' 
#' @param x The single-cell data object.
#' @param assayType The name of the assay slot to access. Not required for 
#'   Seurat v1 or v2.
#' @name getExpr
#' @export
#' 
setGeneric("getExpr",function(x,assayType) standardGeneric("getExpr"))


# ^ getMD ----

#' Get metadata from input data object
#'
#' Extract the cell metadata data frame from a single-cell data object
#' containing the input data for scClustViz visualization.
#'
#' This is a wrapper function to the relevant class's cell metadata slot
#' accessor / assignment method. Currently supported input object classes:
#' \itemize{
#'   \item Class \code{\link[Seurat]{seurat}/\link[Seurat]{Seurat}} accessed by 
#'     \code{x@data.info} or \code{x@meta.data}, 
#'     depending on Seurat object version.
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
#'   \item Class \code{\link[Seurat]{seurat}/\link[Seurat]{Seurat}} accessed by 
#'     \code{x@DRtype.rot} or \code{x@dr$DRtype@cell.embeddings} or 
#'     \code{x@reductions$DRtype@cell.embeddings}, 
#'     depending on Seurat object version.
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

# ^ seurat (v1/2) ----
suppressMessages(
  setMethod("getExpr","seurat",function(x) {
    slot(x,"data")
  })
)

suppressMessages(
  setMethod("getMD","seurat",function(x) {
    if (.hasSlot(x,"meta.data")) {
      slot(x,"meta.data")
    } else {
      slot(x,"data.info") #Seurat v1
    }
  })
)

suppressMessages(
  setMethod("getEmb","seurat",function(x,DRtype) {
    if (.hasSlot(x,"dr")) {
      if (missing(DRtype)) {
        stop(paste(paste0("DRtype must be specified."),
                   "The following cell embeddings are available in this object:",
                   paste0(names(slot(x,"dr")),collapse=", "),sep="\n  "))
      }
      if (tolower(DRtype) %in% names(slot(x,"dr"))) {
        slot(x@dr[[tolower(DRtype)]],"cell.embeddings")
      } else {
        stop(paste(paste0("DRtype '",DRtype,"' not found."),
                   "The following cell embeddings are available in this object:",
                   paste0(names(slot(x,"dr")),collapse=", "),sep="\n  "))
      }
    } else {
      if (.hasSlot(x,paste0(tolower(DRtype),".rot"))) {
        slot(x,paste0(tolower(DRtype),".rot"))
      } else {
        oldSrots <- c("pca","ica","tsne")
        oldSrots <- oldSrots[sapply(oldSrots,function(X) .hasSlot(x,paste0(X,".rot")))]
        stop(paste(paste0("DRtype '",DRtype,"' not found."),
                   "The following cell embeddings are available in this object:",
                   paste0(oldSrots,collapse=", "),sep="\n  "))
      }
    }
  })
)



# ^ Seurat (v3) ----

suppressMessages(
  setMethod("getExpr","Seurat",function(x,assayType) {
    if (missing(assayType)) {
      stop(paste(paste0("assayType must be specified."),
                 "The following assay data are available in this object:",
                 paste0(names(slot(x,"assays")),collapse=", "),sep="\n  "))
    }
    if (assayType %in% names(slot(x,"assays"))) {
      return(x@assays[[assayType]]@data)
    } else {
      stop(paste(paste0("assayType '",assayType,"' not found."),
                 "The following assay data are available in this object:",
                 paste0(names(slot(x,"assays")),collapse=", "),sep="\n  "))
    }
  })
)

suppressMessages(
  setMethod("getMD","Seurat",function(x) {
    return(slot(x,"meta.data"))
  })
)

suppressMessages(
  setMethod("getEmb","Seurat",function(x,DRtype) {
    if (missing(DRtype)) {
      stop(paste(paste0("DRtype must be specified."),
                 "The following cell embeddings are available in this object:",
                 paste0(names(slot(x,"reductions")),collapse=", "),sep="\n  "))
    }
    if (tolower(DRtype) %in% names(slot(x,"reductions"))) {
      return(slot(x@reductions[[tolower(DRtype)]],"cell.embeddings"))
    } else {
      stop(paste(paste0("DRtype '",DRtype,"' not found."),
                 "The following cell embeddings are available in this object:",
                 paste0(names(slot(x,"reductions")),collapse=", "),sep="\n  "))
    }
  })
)


# ^ SingleCellExperiment ----
suppressMessages(
  setMethod("getExpr","SingleCellExperiment",function(x,assayType) {
    if (missing(assayType)) {
      stop(paste(paste0("assayType must be specified."),
                 "The following assay data are available in this object:",
                 paste0(SummarizedExperiment::assayNames(x),collapse=", "),
                 sep="\n  "))
    }
    if (assayType %in% SummarizedExperiment::assayNames(x)) {
      return(SummarizedExperiment::assay(x,assayType))
    } else {
      stop(paste(paste0("assayType '",assayType,"' not found."),
                 "The following assay data are available in this object:",
                 paste0(SummarizedExperiment::assayNames(x),collapse=", "),
                 sep="\n  "))
    }
  })
)

suppressMessages(
  setMethod("getMD","SingleCellExperiment",
            function(x) SingleCellExperiment::colData(x))
)

suppressMessages(
  setMethod("getEmb","SingleCellExperiment",function(x,DRtype) { 
    if (missing(DRtype)) {
      stop(paste(paste0("DRtype must be specified."),
                 "The following cell embeddings are available in this object:",
                 paste0(SingleCellExperiment::reducedDimNames(x),collapse=", "),
                 sep="\n  "))
    }
    if (tolower(DRtype) %in% tolower(SingleCellExperiment::reducedDimNames(x))) {
      return(SingleCellExperiment::reducedDim(
        x,
        SingleCellExperiment::reducedDimNames(x)[
          tolower(DRtype) == tolower(SingleCellExperiment::reducedDimNames(x))
          ]))
    } else {
      stop(paste(paste0("DRtype '",DRtype,"' not found."),
                 "The following cell embeddings are available in this object:",
                 paste0(SingleCellExperiment::reducedDimNames(x),collapse=", "),
                 sep="\n  "))
    }
  })
)