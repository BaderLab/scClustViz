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
#' accessor method. Currently supported input object classes: \itemize{ \item
#' Class \code{\link[Seurat]{seurat}/\link[Seurat]{Seurat}} stored in
#' \code{x@data} or \code{x@assays[[assayType]]@assaySlot}, depending on Seurat
#' object version. \item Class
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} accessed by
#' \code{\link[SummarizedExperiment]{assay}(x,assayType)}. }
#' \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#' for other data objects here!}
#'
#' @param x The single-cell data object.
#' @param assayType A length-one character vector representing the assay object
#'   in which the expression data is stored in the input object. For Seurat v1
#'   or v2 objects, set this to "". For Seurat v3 objects, this is often "RNA".
#'   For SingleCellExperiment objects, this is often "logcounts". See Details
#'   for how this argument is used in the accessor functions for each class.
#' @param assaySlot An optional length-one character vector representing the
#'   slot of the Seurat v3 \code{\link[Seurat]{Assay}} object to use. In Seurat
#'   v3, normalized data is stored in the "data" slot, and counts in the
#'   "counts" slot. See Details for how this argument is used in the accessor
#'   functions for each class.
#' @name getExpr
#' @export
#' 
setGeneric("getExpr",function(x,assayType,assaySlot) standardGeneric("getExpr"))


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


# ^ hasEmb ----

#' List cell embedding types from input data object
#'
#' Show all available cell embedding coordinates from the dimensionality reduction
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
#' @name hasEmb
#' @export
#' 
setGeneric("hasEmb",function(x) standardGeneric("hasEmb"))


# ^ subsetCells ----

#' Subset a data object to return just the specified cells
#' 
#' Internal fx.  Way shittier than the modern object class' subset fx. 
#' Req'd because old Seurat is a POS.
#' 
#' This is a wrapper function to the relevant class's subsetting method. 
#' Currently supported input object classes:
#' \itemize{
#'   \item Class \code{\link[Seurat]{seurat}/\link[Seurat]{Seurat}}.
#'   \item Class \code{\link[SingleCellExperiment]{SingleCellExperiment}}.
#' }
#' \href{https://github.com/BaderLab/scClustViz/issues}{Please submit requests
#' for other data objects here!}
#'
#' @param x The single-cell data object.
#' @param cells A character vector of cell names to extract. This is 
#'   intentionally less versatile than most subsetting, since older versions of 
#'   Seurat required cell names explicitly.
#' @name subsetCells
#' 
setGeneric("subsetCells",function(x,cells) standardGeneric("subsetCells"))


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
      if (DRtype %in% names(slot(x,"dr"))) {
        slot(x@dr[[DRtype]],"cell.embeddings")
      } else {
        stop(paste(paste0("DRtype '",DRtype,"' not found."),
                   "The following cell embeddings are available in this object:",
                   paste0(names(slot(x,"dr")),collapse=", "),sep="\n  "))
      }
    } else {
      if (.hasSlot(x,paste0(DRtype,".rot"))) {
        slot(x,paste0(DRtype,".rot"))
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


suppressMessages(
  setMethod("hasEmb","seurat",function(x) {
    if (.hasSlot(x,"dr")) {
      return(names(slot(x,"dr")))
    } else {
      oldSrots <- c("pca","ica","tsne")
      oldSrots <- oldSrots[sapply(oldSrots,function(X) .hasSlot(x,paste0(X,".rot")))]
      return(oldSrots)
    }
  })
)


suppressMessages(
  setMethod("subsetCells","seurat",function(x,cells) {
    Seurat::SubsetData(x,cells.use=cells)
  })
)


# ^ Seurat (v3) ----
suppressMessages(
  setMethod("getExpr","Seurat",function(x,assayType,assaySlot) {
    if (missing(assayType)) {
      stop(paste(paste0("assayType must be specified."),
                 "The following assay data are available in this object:",
                 paste0(names(slot(x,"assays")),collapse=", "),sep="\n  "))
    }
    if (assayType %in% names(slot(x,"assays"))) {
      if (!missing(assaySlot)) {
        if (is.na(assaySlot) | assaySlot == "") {
          return(x@assays[[assayType]]@data)
        } else {
          return(slot(x@assays[[assayType]],assaySlot))
        }
      } else {
        return(x@assays[[assayType]]@data)
      }
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
    if (DRtype %in% names(slot(x,"reductions"))) {
      return(slot(x@reductions[[DRtype]],"cell.embeddings"))
    } else {
      stop(paste(paste0("DRtype '",DRtype,"' not found."),
                 "The following cell embeddings are available in this object:",
                 paste0(names(slot(x,"reductions")),collapse=", "),sep="\n  "))
    }
  })
)


suppressMessages(
  setMethod("hasEmb","Seurat",function(x) {
    return(names(slot(x,"reductions")))
  })
)


suppressMessages(
  setMethod("subsetCells","Seurat",function(x,cells) {
    subset(x,cells=cells)
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
    if (DRtype %in% SingleCellExperiment::reducedDimNames(x)) {
      return(SingleCellExperiment::reducedDim(
        x,
        SingleCellExperiment::reducedDimNames(x)[
          DRtype == SingleCellExperiment::reducedDimNames(x)
          ]))
    } else {
      stop(paste(paste0("DRtype '",DRtype,"' not found."),
                 "The following cell embeddings are available in this object:",
                 paste0(SingleCellExperiment::reducedDimNames(x),collapse=", "),
                 sep="\n  "))
    }
  })
)


suppressMessages(
  setMethod("hasEmb","SingleCellExperiment",function(x) { 
    SingleCellExperiment::reducedDimNames(x)
  })
)


suppressMessages(
  setMethod("subsetCells","SingleCellExperiment",function(x,cells) {
    return(x[,cells])
  })
)
